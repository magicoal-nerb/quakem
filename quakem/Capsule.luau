--!strict

-- handles general roblox capsule collision physics
-- deterministic, easy to use for server auth or custom
-- player physics if needed.

-- this uses a custom physics + collision solver in here,
-- where i use position based dynamics for the static capsule collision
-- response and a projected gauss seidel solver for the generic unanchored collisions.

-- note that the collision solver approximates the response because roblox cannot easily
-- get collision data, so i went with speculative contacts for CCD (which is guaranteed
-- to not clip through stuff), and it handles deep contacts aswell (through using an
-- approximate box SAT and checking outside perimeter.) however, in practice, these
-- should be unnoticeable from using the proven algs (gjk/sat/mpr/dcol.) so, if you
-- already have those algs suited for your assets, then you can replace mine with those.

-- this also supports generic mesh part collision (because roblox's get closest point method
-- supports it), and it supports terrain (which internally uses a mini GJK implementation
-- for getting closest points of arbritrary terrain vertices)

-- for generic mesh part collision, beware of concave geometry because otherwise
-- normal reconstruction can fail (a box strip for example, may not recognize a wall properly
-- because of the nature of the solver using GetClosestPointOnSurface. if this happens,
-- just disable collision of that meshpart and reconstruct its geometry with primitives instead.)

-- magicoal_nerb :3

local Terrain = workspace:WaitForChild("Terrain")
local Capsule = {}

export type Feature = { delta: Vector3, pA: Vector3, pB: Vector3 }
export type Contact = {
	normal: Vector3,
	point: Vector3,
	part: BasePart,
	mtv: Vector3,
}

export type Capsule = {
	position: Vector3,
	velocity: Vector3,
	angularVelocity: Vector3,

	height: number,
	radius: number,
	friction: number,
	
	mass: number,
	size: Vector3,
	
	contacts: { Contact },
	rotation: CFrame,
	frame: CFrame,

	overlap: OverlapParams,
	params: RaycastParams,
}

-- How similar should a normal be before we prune it
const kPlaneMerge = 0.9

-- Friction
const kGroundFriction = 14

-- Bias for contact depth
const kSolverBias = 30

-- Arbritrary occupancy cutoff for considering air/ground
const kTerrainCutoff = 0.2

-- Number of iterations per terrain GJK call
const kGjkIters = 32

const kTerrainVerticies: { Vector3 } = {
	Vector3.new(0, 1, 0), Vector3.new(1, 1, 0), Vector3.new(0, 1, 1), Vector3.new(1, 1, 1),
	Vector3.new(0, 0, 0), Vector3.new(1, 0, 0), Vector3.new(0, 0, 1), Vector3.new(1, 0, 1)
}

const kTerrainEdges: { number } = {
	0x15, 0x26, 0x37, 0x48,
	0x56, 0x78, 0x57, 0x68,
	0x12, 0x34, 0x13, 0x24
}

local function utilInverseInertia(R: CFrame, size: Vector3, mass: number): CFrame
	-- Calculates the inverse moment of inertia tensor of a cuboid. We assume
	-- this for all part types, since it is a good estimate. Roblox does not
	-- provide us information about the local inertia tensor.

	-- L_world = RL_local => R^T L_world ; w_world = Rw_local => R^T w_world
	-- <=> R^T L_world = I_local R^T w_world
	-- <=> L_world = R I_local R^T w_world (since R^-1 = R^T)
	-- So, I = R I_local R^T.
	-- I^-1 = (R^T)^-1 * I_local^-1 * R^-1 <=> R * I_local^-1 * R^T

	local sx, sy, sz = size.X, size.Y, size.Z
	local localInertia = Vector3.new(sy*sy + sz*sz, sx*sx + sz*sz, sx*sx + sy*sy) * (mass / 12)
	local invLocalInertia = CFrame.new(0, 0, 0, 1/localInertia.X, 0, 0, 0, 1/localInertia.Y, 0, 0, 0, 1/localInertia.Z)
	return R * invLocalInertia * R:Inverse()
end

local function utilConstrain(disp: Vector3, normal: Vector3): Vector3
	-- Minimize dot(disp, normal), given velocity + normal*lambda_1
	local dot = math.min(disp:Dot(normal), 0.0)
	return disp - normal * dot
end

local function utilRelevantPlanes(capsule: Capsule, planes: { Contact }): { Contact }
	-- This approximates an actual active-set solver solution to N^T lambda >= 0,
	-- given that lambda >= 0. However, when I implemented the technique, the technique
	-- sometimes overcorrected compared to the greedy approach.
	table.sort(planes, function(a: Contact, b: Contact)
		return a.normal:Dot(a.mtv) > b.normal:Dot(b.mtv)
	end)

	-- Remove duplicate entries (we're choosing the most impactful collisions
	-- to the least impactful greedily.)
	local active = {} :: { Contact }
	for _, plane in planes do
		-- Discard planes where we are moving away from them or if the part
		-- associated with them is unanchored
		if plane.normal:Dot(capsule.velocity) > 0.0 then continue end
		if plane.normal:Dot(plane.mtv) <= 0.0 then continue end
		if plane.part == Terrain then continue end
		if not plane.part.Anchored then continue end
		
		local duplicate = false
		for _, other in active do
			if plane.normal:Dot(other.normal) > kPlaneMerge then
				duplicate = true
				break
			end
		end

		if not duplicate and #active < 3 then
			table.insert(active, plane)
		end
	end

	return active	
end

local function utilTerrainCell(sampled: { number }): { Vector3 }
	-- Gets the terrain cell vertices. This returns a set of vertices
	-- that are on this cell.
	local vertices = {} :: { Vector3 }
	
	for _, edge in kTerrainEdges do
		-- One of the sides must be occupied with air
		-- if we want to compute the mesh.
		local va = bit32.band(edge, 0xF)
		local vb = bit32.rshift(edge, 4)

		local a = sampled[va]
		local b = sampled[vb]
		if a > b then
			-- Need a < b.
			va, vb = vb, va;
			a, b = b, a;
		end

		if a < kTerrainCutoff and b >= kTerrainCutoff then
			-- There's an intersection at this edge, so we
			-- add it to the convex hull
			local t = (kTerrainCutoff - a) / (b - a)
			local p = kTerrainVerticies[va]:Lerp(kTerrainVerticies[vb], t) + Vector3.one * 0.5
			table.insert(vertices, 4 * p)
		end
	end

	return vertices
end

local function utilLinePoint(a: Vector3, b: Vector3, p: Vector3): (Vector3, Vector3)
	-- dot(a + ab*t, ab) = dot(p, ab)
	local ab = b - a
	local abab = ab:Dot(ab)

	if abab < 1e-6 then return a, Vector3.xAxis end

	local t = math.clamp((p - a):Dot(ab) / abab, 0.0, 1.0)
	return a + ab * t, Vector3.new(1.0 - t, t, 0)
end

local function utilClosestPointTri(a: Vector3, b: Vector3, c: Vector3, p: Vector3): (Vector3, Vector3)
	-- Taken from Christer Ericson's Real Time Collision book
	local ab, ac, ap, bp, cp = b - a, c - a, p - a, p - b, p - c
	local norm = ab:Cross(ac)

	if norm:Dot(norm) < 1e-3 then
		-- If it's degenerate, then check the possible lines
		local ab, wAB = utilLinePoint(a, b, p) ab -= p
		local ac, wAC = utilLinePoint(a, c, p) ac -= p
		local bc, wBC = utilLinePoint(b, c, p) bc -= p

		local abab, acac, bcbc = ab:Dot(ab), ac:Dot(ac), bc:Dot(bc)
		if abab < acac and abab < bcbc then
			return ab + p, wAB
		elseif acac < bcbc then
			return ac + p, Vector3.new(wAC.X, 0, wAC.Y)
		else
			return bc + p, Vector3.new(0, wBC.X, wBC.Y)
		end
	end

	local d1, d2, d3, d4, d5, d6 = ab:Dot(ap), ac:Dot(ap), ab:Dot(bp), ac:Dot(bp),
		ab:Dot(cp), ac:Dot(cp)

	if d1 <= 0.0 and d2 <= 0.0 then
		-- Inside of AB and AC half plane, so choose A 
		return a, Vector3.new(1, 0, 0)
	elseif d3 >= 0.0 and d4 <= d3 then
		-- Then, we are on the half plane facing AB and BP, so
		-- we choose B
		return b, Vector3.new(0, 1, 0)
	elseif d6 >= 0.0 and d5 <= d6 then
		-- Then we are on the half plane facing AC and CP, so we
		-- choose C
		return c, Vector3.new(0, 0, 1)
	end

	-- Calculate the determinant of the D entries, which is actually just a cross
	-- product in disguise (see their derivation in the book for more info)
	local vc, vb, va = d1 * d4 - d3 * d2, d5 * d2 - d1 * d6, d3 * d6 - d5 * d4
	if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0 then
		-- Then it's just AB
		local v = d1 / (d1 - d3)
		return a + v * ab, Vector3.new(1.0 - v, v, 0)
	elseif vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0 then
		-- Then it's just AC
		local w = d2 / (d2 - d6)
		return a + w * ac, Vector3.new(1.0 - w, 0, w)
	elseif va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0 then
		-- Then it's just BC
		local w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
		return b + w * (c - b), Vector3.new(0, 1.0 - w, w)
	end

	-- Interior of poly
	local denom = 1.0 / (va + vb + vc)
	local v, w = vb * denom, vc * denom
	return a + ab * v + ac * w, Vector3.new(1.0 - v - w, v, w)
end

local function utilGjkClosestPoint(suppA: Vector3, suppB: (Vector3) -> Vector3, delta: Vector3): Vector3
	local support = function(dir: Vector3): Feature
		local pA, pB = suppA, suppB(dir)
		return {
			delta = pA - pB,
			pA = pA, pB = pB,
		}
	end

	-- Hopefully a pretty neat approach for doing GJK, where we just continually do the tetrahedron
	-- case. Seems to converge similarly to the standard approach (see Capsule.luau for that.)
	local a = support(delta)
	local b = support(-a.delta)
	local c = support(-utilLinePoint(a.delta, b.delta, Vector3.zero))

	for i = 1, kGjkIters do
		local point = utilClosestPointTri(a.delta, b.delta, c.delta, Vector3.zero)
		if point:FuzzyEq(Vector3.zero) then
			-- We hit the origin, so just return the zero vector
			break
		end

		-- Search from the closest point for the new point
		local normal = point.Unit
		local d = support(normal)

		local dot = d.delta:Dot(normal)
		local wA, wB, wC = a.delta:Dot(normal), b.delta:Dot(normal), c.delta:Dot(normal)
		if dot >= math.min(wA, wB, wC) then
			-- Converged because our distance from the origin did not improve
			-- (GJK closest point distance has to decrease monotonically.)
			break
		end

		-- Replace the furthest triangle (this is basically the tetrahedron case
		-- but we try to decrease our distance the most!)
		local wACD = utilClosestPointTri(a.delta, c.delta, d.delta, Vector3.zero)
		local wABD = utilClosestPointTri(a.delta, b.delta, d.delta, Vector3.zero)
		local wBCD = utilClosestPointTri(b.delta, c.delta, d.delta, Vector3.zero)

		local acd = wACD:Dot(wACD)
		local abd = wABD:Dot(wABD)
		local bcd = wBCD:Dot(wBCD)
		local abc = point:Dot(point)

		if acd <= abd and acd <= bcd and acd <= abc then
			-- Choose ACD
			b = d
		elseif abd <= bcd and abd <= abc then
			-- Choose ABD
			c = d
		elseif bcd <= abc then
			-- Choose BCD
			a = d
		else
			break
		end
	end

	local point, barycoord = utilClosestPointTri(a.delta, b.delta, c.delta, Vector3.zero)
	return a.pB * barycoord.X + b.pB * barycoord.Y + c.pB * barycoord.Z
end

local function utilSolveActive(capsule: Capsule, planes: { Contact })
	if #planes == 0 then
		-- Nothing
		return
	elseif #planes == 1 then
		local a = planes[1]

		-- Push away from nA
		local nA = a.normal
		local d0, v0 = nA:Dot(a.mtv), nA:Dot(capsule.velocity)

		capsule.velocity -= nA * v0
		capsule.position += nA * d0
	elseif #planes == 2 then
		local a, b = planes[1], planes[2]

		-- The only real direction we can move here is orthogonal to nA and nB
		local nA, nB = a.normal, b.normal
		local nAxB = nA:Cross(nB)

		local det = nAxB:Dot(nAxB)
		if math.abs(det) < 1e-6 then return end

		local d0, v0 = nA:Dot(a.mtv), nA:Dot(capsule.velocity)
		local d1, v1 = nB:Dot(b.mtv), nB:Dot(capsule.velocity)

		capsule.position += (nB:Cross(nAxB) * d0 + nAxB:Cross(nA) * d1) / det
		capsule.velocity -= (nB:Cross(nAxB) * v0 + nAxB:Cross(nA) * v1) / det
	else
		-- Concatenate all of the normal vectors
		-- as a column vector and solve by inverting the matrix 
		local a, b, c = planes[1], planes[2], planes[3]
		local nA, nB, nC = a.normal, b.normal, c.normal

		local det = nC:Dot(nA:Cross(nB))
		if math.abs(det) < 1e-6 then return end

		local d0, d1, d2 = nA:Dot(a.mtv), nB:Dot(b.mtv), nC:Dot(c.mtv)
		local v0, v1, v2 = nA:Dot(capsule.velocity), nB:Dot(capsule.velocity), nC:Dot(capsule.velocity)
		capsule.position += (d0 * nB:Cross(nC) + d1 * nC:Cross(nA) + d2 * nA:Cross(nB)) / det
		capsule.velocity -= (v0 * nB:Cross(nC) + v1 * nC:Cross(nA) + v2 * nA:Cross(nB)) / det
	end
end

function Capsule.push(collider: Capsule, what: BasePart, mtv: Vector3, linFac: number)
	if mtv:FuzzyEq(Vector3.zero) then return end

	-- We assume a typical PGS formulation.
	-- Also, this is exactly how you'd do a contact constraint
	-- in a typical physics engine :3
	local pb = what:GetClosestPointOnSurface(collider.position)
	local rA, rB = pb - collider.position, pb - what.Position

	-- Calculate the friction directions
	local normal = mtv.Unit
	local up = math.abs(normal.Y) < 0.99 and Vector3.yAxis or Vector3.xAxis
	local tangent = up:Cross(normal).Unit
	local bitangent = normal:Cross(tangent).Unit

	-- Get the current state variables
	local vB, wB = what.AssemblyLinearVelocity, what.AssemblyAngularVelocity
	local vA, wA = collider.velocity, collider.angularVelocity

	local cframeA, cframeB = CFrame.new(collider.position) * collider.rotation, what.CFrame
	local sA, sB = Vector3.new(collider.radius, collider.height, collider.radius) * 2, what.Size
	local iA, iB = utilInverseInertia(collider.rotation, sA, collider.mass), utilInverseInertia(cframeB.Rotation, sB, what.AssemblyMass)
	local mA, mB = 1.0 / collider.mass, 1.0 / what.AssemblyMass

	-- Precompute cross products
	local rnA, rtA, rbA = rA:Cross(normal), rA:Cross(tangent), rA:Cross(bitangent)
	local rnB, rtB, rbB = rB:Cross(normal), rB:Cross(tangent), rB:Cross(bitangent)

	-- Calculate the effective mass, which is from J M^-1 J^T
	-- see: https://github.com/magicoal-nerb/impulse/blob/main/luau-physics/narrowphase/Contact.luau
	local effMass = mA + mB + rnB:Dot(iB * rnB) + rnA:Dot(iA * rnA)
	local effMassT = mA + mB + rtB:Dot(iB * rtB) + rtA:Dot(iA * rtA)
	local effMassB = mA + mB + rbB:Dot(iB * rbB) + rbA:Dot(iA * rbA)

	if effMass < 1e-8 or effMassT < 1e-8 or effMassB < 1e-8 then
		-- This would go terribly wrong
		return
	end

	-- Relative velocity is rA's velocity - rB's velocity.
	-- Enforce dot(relativeVelocity, currentSurface) >= 0.0
	local invEffMass = 1.0 / effMass
	local rV = vA + wA:Cross(rA) - vB - wB:Cross(rB)

	local nDot = rV:Dot(normal) - math.max(mtv.Magnitude * kSolverBias, 0.0)
	local tDot, bDot = -rV:Dot(tangent), -rV:Dot(bitangent)
	local maxFriction = math.max(-nDot * invEffMass, 0.0) * collider.friction

	-- Make a simplistic friction square (since, we don't need to do anything fancy
	-- for a rotation effect, lol.) Ideally, you'd use a 2x2 matrix inverse for this.
	local tImpulse = math.clamp(tDot / effMassT, -maxFriction, maxFriction)
	local bImpulse = math.clamp(bDot / effMassB, -maxFriction, maxFriction)
	local friction = tangent * tImpulse + bitangent * bImpulse

	-- Apply impulse on collider and candidate. Do not apply friction
	-- onto the collider's linear velocity (but angular velocity is okay.)
	local impulse = math.max(-nDot * invEffMass, 0.0) * normal + friction
	collider.angularVelocity += iA * rA:Cross(impulse)

	collider.position += mtv * linFac
	collider.velocity += mA * (impulse - friction) * linFac

	if not what.Anchored then
		what.AssemblyLinearVelocity -= impulse * mB
		what.AssemblyAngularVelocity -= iB * rB:Cross(impulse)
	end
end

function Capsule.getContacts(capsule: Capsule, change: Vector3): { BasePart }
	-- Gets all possible parts that could contact with our
	-- given capsule.
	local height, radius = capsule.height, capsule.radius
	local halfSize = Vector3.new(radius, height, radius)
	
	local searchSize = (halfSize + change:Abs()) * 2	
	local results = workspace:GetPartBoundsInBox(
		CFrame.new(capsule.position + change),
		searchSize,
		capsule.overlap
	)
	
	-- Appends fake BaseParts in `parts` such that we can handle
	-- collisions with the terrain data.
	local initPos = capsule.position + change
	local reg = Region3.new(
		initPos - searchSize * 0.5 - Vector3.one * 4,
		initPos + searchSize * 0.5 + Vector3.one * 4
	):ExpandToGrid(4)

	debug.profilebegin("quakec::get-terrain-contacts")

	-- Sample adjacent cells for clean voxel collision!
	local mats, occs = Terrain:ReadVoxels(reg, 4)
	local size = mats.Size
	local contacts = table.create(8, 0.0)

	local sample = function(vec: Vector3)
		if not (vec - size):Max(Vector3.zero):FuzzyEq(Vector3.zero) then
			-- Out of bounds!
			return 0.0
		end

		return occs[vec.X + 1][vec.Y + 1][vec.Z + 1]
	end

	size -= Vector3.one
	for i = 0, size.X * size.Y * size.Z - 1 do
		-- Generate the terrain cell dynamically
		-- from the lookup table.
		local vec = Vector3.new(i // (size.Y * size.X), (i // size.Z) % size.Y, i % size.Z)
		for j = 1, #kTerrainVerticies do
			contacts[j] = sample(vec + kTerrainVerticies[j])
		end

		local centroid = reg.CFrame.Position - reg.Size * 0.5 + vec * 4
		local result = utilTerrainCell(contacts)
		if #result == 0 then
			-- It's air, so don't bother.
			continue
		end

		-- Create a support mapping of this object
		local support = function(dir: Vector3)
			local bestDot = -math.huge
			local bestIdx = 1

			for j, point in result do
				local dot = point:Dot(dir)
				if dot > bestDot then
					bestDot = dot
					bestIdx = j
				end
			end

			return result[bestIdx] + centroid
		end

		-- This could be improved to not create a table
		-- through caching terrain cells. But for now, we'll
		-- just take the temporary memory penalty for the sake
		-- of having a clean solver.
		table.insert(results, {
			CFrame = CFrame.new(centroid),
			Size = Vector3.one * 4,
			Anchored = true,
			AssemblyMass = math.huge,
			AssemblyLinearVelocity = Vector3.zero,
			AssemblyAngularVelocity = Vector3.zero,
			GetClosestPointOnSurface = function(_, pos: Vector3)
				return utilGjkClosestPoint(pos, support, pos - centroid)
			end,
		} :: BasePart)
	end

	debug.profileend()
	
	return results
end

function Capsule.collide(c: Capsule, dt: number)
	-- Handles collision logic
	local size = Vector3.new(c.radius, c.height, c.radius) * 2	
	local height, radius = c.height, c.radius
	height -= radius
	
	debug.profilebegin("capsule::speculative-solver")

	-- Check if the inner portion of the capsule
	-- collided with anything. Note that in Cylinder, I actually
	-- solve this portion analytically, and uses a parabolic trace.
	local result = workspace:Blockcast(
		CFrame.new(c.position),
		size - Vector3.one * radius * 2,
		c.velocity * dt,
		c.params
	)
	
	if result then
		-- We found a contact! Simply just resolve it.
		local newPosition, normal = result.Position, result.Normal
		local point = newPosition - normal * radius
		local delta = c.velocity * dt
		
		c.position = result.Position
		c.velocity = utilConstrain(c.velocity, normal)
		dt *= 1.0 - math.clamp((newPosition - c.position):Dot(delta) / delta:Dot(delta), 0, 1)
		Capsule.push(c, result.Instance :: BasePart, Vector3.zero, 0.0)
	end
	

	local candidates = Capsule.getContacts(c, c.velocity * dt)
	c.position += c.velocity * dt

	debug.profileend()

	debug.profilebegin("capsule::push-solver")

	for i = 1, 2 do
		table.clear(c.contacts)
		
		-- We alternate the heights so it works on concave shapes
		local offset = Vector3.new(0, math.lerp(-height, height, i % 2), 0)
		
		for j, candidate in candidates do
			-- We estimate the collision candidate through calculating a possible
			-- plane of intersection.
			local pb = candidate:GetClosestPointOnSurface(c.position + offset)
			local pa = utilLinePoint(c.position + Vector3.new(0, height, 0), c.position - Vector3.new(0, height, 0), pb)

			local mtv = pa - pb
			local diff = radius - mtv.Magnitude
			
			if diff < -1e-4 then
				-- Not colliding. Also this slop has to exist
				-- because it causes MeshParts to hover (although, primitives work fine.)
				continue
			elseif mtv:Dot(mtv) < 1e-4 then
				-- Capsule is completely engulfed
				local cframe, size = candidate.CFrame, candidate.Size

				local absCFrame = CFrame.fromMatrix(cframe.Position, cframe.XVector:Abs(), cframe.YVector:Abs(), cframe.ZVector:Abs())
				local min, max = absCFrame * -size, absCFrame * size

				-- We'll just do a crude SAT box check and move them out of it
				-- (since we can't really assume much from general geometry.)
				local push = (c.position - min):Min(max - c.position):Max(Vector3.zero)
				local normal = Vector3.zero
				
				if push.X < push.Y and push.X < push.Z then
					normal = Vector3.xAxis * math.sign(c.position.X - cframe.Position.X)
				elseif push.Y < push.Z then
					normal = Vector3.yAxis * math.sign(c.position.Y - cframe.Position.Y)
				else
					normal = Vector3.zAxis * math.sign(c.position.Z - cframe.Position.Z)
				end

				mtv = candidate:GetClosestPointOnSurface(c.position + normal * push) - c.position
				mtv += normal * normal:Dot(size * normal:Sign() * 0.5)
				
				normal = mtv.Unit
				diff = mtv.Magnitude
			end
			
			local normal: Vector3 = mtv.Unit
			if typeof(candidate) == "Instance" then
				Capsule.push(c, candidate, diff * normal, candidate.Anchored and 0.0 or 1.0)
			else
				c.position += diff * normal
				c.velocity = utilConstrain(c.velocity, normal)
			end
			
			-- Append this contact here!
			table.insert(c.contacts, {
				mtv = diff * normal,
				normal = normal,
				point = pb,
				part = typeof(candidate) == "Instance" and candidate or Terrain,
			} :: Contact)
		end
		
		utilSolveActive(c, utilRelevantPlanes(c, c.contacts))
	end

	debug.profileend()
end

Capsule.utilInverseInertia = utilInverseInertia

return Capsule