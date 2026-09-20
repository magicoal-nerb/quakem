--!strict

-- cylinder cylinder, you can also easily make this into a capsule cylinder
-- by changing the utilCylinderSupport function.

-- most notable parts are the physics solver + collision response. physics solver is a traditional
-- PGS solver, but one sided against the physics object being interacted with.
-- the continuous collision part uses an analytical TOI for resolution (respects velocity + acceleration,
-- the parabolic path.). contact resolution is an analytical approximation of the active-set approach through
-- greedily selecting contact depth and normal dedup, which helps satisfy the MLCP
-- N^T x >= 0, x >= 0. using this approach seems to fix the active-set overcorrecting the solution

-- main limitation is that this mainly works on convex geometry, so concave mesh parts
-- will not work as cleanly. also, this is deterministic and analytically solves
-- all of the primitives, where you get an exact TOI from a given velocity and acceleration.

-- uses a probably novel idea where we just use the closest points for the separating axes
-- for the cylinder (mostly just because we exploit the rotational symmetry from a cylinder/capsule shape.)
-- moreover, the SAT is modified to also support continuous collision, which should be faster than FEV/VClip
-- analytical methods for basic primitives.

-- this supports solving any general SAT configuration, and general convex
-- support functions for GJK. either of which work cleanly hopefully too :3
-- magicoal_nerb <3

export type Contact = { position: Vector3, normal: Vector3, mtv: Vector3, part: BasePart }
const kPlaneMerge = 0.9
const kSolverBias = 30

export type Feature = { delta: Vector3, pA: Vector3, pB: Vector3 }
export type Support = (Vector3) -> (Vector3)
export type Cylinder = {
	position: Vector3,
	angularVelocity: Vector3,
	acceleration: Vector3,
	velocity: Vector3,

	rotation: CFrame,
	mass: number,
	friction: number,

	height: number,
	radius: number,
	params: RaycastParams,
	overlap: OverlapParams,
	contacts: { Contact },
}

const kGjkIters = 4
const kSatIters = 4

local Cylinder = {}

local function utilRayPlane(p: Vector3, v: Vector3, s: Vector3, n: Vector3): (Vector3, number)
	-- dot(p + t * v, n) = dot(s, n)
	-- <=> t = dot(s - p, n) / dot(v, n)
	local ps = s - p
	local t = ps:Dot(n) / v:Dot(n)
	return p + t * v, t
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

	-- Calculate the determinant of [d1 d2 ; d3 d4] 
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

local function utilGjkClosestPoint(suppA: (Vector3) -> Vector3, suppB: (Vector3) -> Vector3, delta: Vector3): (Vector3, Vector3)
	local support = function(dir: Vector3): Feature
		local pA, pB = suppA(-dir), suppB(dir)
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
	return a.pA * barycoord.X + b.pA * barycoord.Y + c.pA * barycoord.Z,
		a.pB * barycoord.X + b.pB * barycoord.Y + c.pB * barycoord.Z
end

local function utilCylinderSupport(cylinder: Cylinder): (Vector3) -> Vector3
	-- Creates a support function for the cylinder I guess
	local radius, height = cylinder.radius, cylinder.height
	return function(dir: Vector3)
		-- Calculate support of the cylinder
		local dt, cc = dir.Y, Vector3.zero
		local pos = cylinder.position
		if math.abs(dt) < 1e-6 then
		    return dir * radius + pos
		else
		    local origin = Vector3.yAxis * height * math.sign(dt)
		    local plane, t = utilRayPlane(Vector3.zero, dir, origin, Vector3.yAxis)
		    local op = plane - origin
		    if op:Dot(op) > radius * radius then
		        op = op.Unit * radius
		    end

		    return origin + op + pos
		end
	end
end

local function utilSolve(a: number, b: number, c: number): number?
	if math.abs(a) < 1e-6 then
		-- Linear case
		if math.abs(b) < 1e-6 then return nil end

		local toi = -c / b
		return toi >= 0.0 and toi or nil
	end

	-- Solve the quadratic
	local disc = b*b - 4*a*c
	if disc < 0.0 then return nil end

	local sqrtDisc = math.sqrt(disc)
	local t0, t1 = (-b - sqrtDisc) / (2*a), (-b + sqrtDisc) / (2*a)

	-- There should only be one minimum positive root
	if t0 >= 0.0 then return t0 end
	if t1 >= 0.0 then return t1 end

	return nil
end

local function utilGenericSat(
	suppA: (Vector3) -> Vector3,
	suppB: (Vector3) -> Vector3,
	velocity: Vector3, accel: Vector3, dt: number,
	sat: { Vector3 }
): (Vector3, number, Vector3, number)
	local push, normal, toiNormal = 1e6, Vector3.zero, Vector3.zero
	local l = 0.0

	for i = 1, kSatIters do
		local pushed = true
		push, normal = 1e6, Vector3.zero

		-- Current state
		local p, v = 0.5 * accel * l * l + velocity * l, velocity + accel * l
		local t = dt - l
		for j, dir in sat do
			-- If dir is NaN
			if dir ~= dir then continue end

			-- Calculate the SAT bounds from A and B
			-- [ min(da,db,dc), max(da,db,dc) ]
			local cMin, cMax = (suppA(-dir) + p):Dot(dir), (suppA(dir) + p):Dot(dir)
			local tMin, tMax = suppB(-dir):Dot(dir), suppB(dir):Dot(dir)

			local d0, d1 = tMin - cMax, cMin - tMax
			if d0 > 0.0 or d1 > 0.0 then
				pushed, push = false, 0.0
				local minD = math.max(d0, d1)
				local dd = v:Dot(dir)
				
				-- Check if capsule/primitive is ahead of eachother
				if d0 > 0 and dd <= 0.0 then continue end
				if d0 < 0 and dd >= 0.0 then continue end

				local ct = utilSolve(
					0.5 * accel:Dot(dir),
					v:Dot(dir),
					d0 > 0.0 and d0 or d1
				)
				
				if not ct then continue end
				if t < ct then continue end

				-- Derivative check, if we're facing the direction
				-- vector, then discard.
				local newVel = v + accel * ct
				local sign = dd < 0.0 and 1 or -1

				if newVel:Dot(dir * sign) >= 0.0 then continue end

				toiNormal = dir * sign
				t = ct
			else
				-- Discrete check
				local minD = math.min(-math.min(d0, 0.0), -math.min(d1, 0.0))
				if minD > push then continue end

				local sign = d0 < d1 and -1 or 1
				normal = -dir * sign
				push = minD
			end
		end

		if pushed then
			-- If we had a discrete contact, we're done!
			break
		elseif t == 0.0 then
			-- We converged on a proper TOI
			break
		end

		l += t
	end

	return normal, push, toiNormal, l
end

local function utilSupport(part: BasePart): (Enum.PartType?, (Vector3) -> Vector3)
	-- Optimization where we just have our own computed support functions
	-- if it is a simple shape. Otherwise, we will just fallback to Roblox's
	-- closest point feature to approximate support functions
	local cframe, halfSize = part.CFrame, part.Size * 0.5
	local shape = part:IsA("Part") and part.Shape or nil
	
	-- Creates a support for this part
	if shape == Enum.PartType.Block then
		return Enum.PartType.Block, function(direction: Vector3)
			local relative = (cframe.Rotation:Inverse() * direction):Sign()
			return cframe:PointToWorldSpace(relative * halfSize)
		end
	elseif shape == Enum.PartType.Wedge or part:IsA("WedgePart") then
		return Enum.PartType.Wedge, function(direction: Vector3)
			local front = Vector3.new(0, halfSize.Z, -halfSize.Y).Unit

			local relative = cframe:VectorToObjectSpace(direction)
			local dot = relative:Dot(front)
			if 0 < dot then
				-- Constrain front face
				local d0 = relative.Y * halfSize.Y + relative.Z * halfSize.Z
				relative = if d0 > 0.0
					then Vector3.new(relative.X, 1, 1)
					else Vector3.new(relative.X, -1, -1)
			end

			return cframe:PointToWorldSpace(relative:Sign() * halfSize)
		end
	elseif shape == Enum.PartType.CornerWedge or part:IsA("CornerWedgePart") then
		return Enum.PartType.CornerWedge, function(direction: Vector3)
			local relative = cframe:VectorToObjectSpace(direction)
			local d = relative * halfSize
			
			-- What's the dot product if we were to check the apex?
			local dotApex = d.X + d.Y - d.Z

			-- Which half plane are we in?
			local bx = if relative.X > 0.0 then 1.0 else -1.0
			local bz = if relative.Z > 0.0 then 1.0 else -1.0

			-- Project dot product onto the best base
			local dotBase = (bx * d.X) - d.Y + (bz * d.Z)
			local finalSign = if dotApex > dotBase
				then Vector3.new(1.0, 1.0, -1.0)
				else Vector3.new(bx, -1.0, bz)

			return cframe:PointToWorldSpace(finalSign * halfSize)
		end
	elseif shape == Enum.PartType.Cylinder then
		local radius = math.min(halfSize.Y, halfSize.Z)
		return Enum.PartType.Cylinder, function(direction: Vector3)
			direction = cframe:VectorToObjectSpace(direction)

			local dotT, cPoint = direction.X, Vector3.new(halfSize.X, 0, 0)
			local h, t, final

			if dotT == 0 then
				-- Just use the inner ring
				final = direction.Unit * radius
			else
				-- Generic case
				cPoint = dotT > 0 and cPoint or -cPoint
				h, t = utilRayPlane(Vector3.zero, direction, cPoint, Vector3.xAxis)
				
				local op = h - cPoint
				if op:Dot(op) > radius * radius then
					op = op.Unit * radius
				end
				
				final = cPoint + op
			end

			return cframe:PointToWorldSpace(final)
		end
	elseif shape == Enum.PartType.Ball then
		local radius = math.min(halfSize.X, halfSize.Y, halfSize.Z)
		return Enum.PartType.Ball, function(direction: Vector3)
			return direction.Unit * radius + cframe.Position
		end
	else
		-- Use a generic support function for this. Will revisit once again
		-- when Roblox adds new shape types
		return nil, function(direction: Vector3)
			return part:GetClosestPointOnSurface(cframe.Position + direction * 1024)
		end
	end
end

local function utilCylinderCastPrimitive(cylinder: Cylinder, part: BasePart, velocity: Vector3, accel: Vector3, dt: number): (Vector3, Vector3, number, Vector3, number)
	-- Generic shape
	local cframe, halfSize = part.CFrame, part.Size * 0.5
	local shape, supportB = utilSupport(part)
	
	local offset = Vector3.yAxis * (cylinder.height - cylinder.radius)
	local pA, pB = utilGjkClosestPoint(supportB, function(dir: Vector3)
		return cylinder.position + math.sign(dir.Y) * offset
	end, part.Position - cylinder.position)
		
	local axes = {
		Vector3.yAxis,
		cframe.XVector, cframe.YVector, cframe.ZVector,
		(pA - pB).Unit,
	}
	
	if shape == Enum.PartType.Wedge then
		-- Add the face normal of the front facing wedge
		local front = Vector3.new(0, halfSize.Z, -halfSize.Y).Unit
		table.insert(axes, cframe.Rotation * front) 
	elseif shape == Enum.PartType.CornerWedge then
		-- Add the face normals of the corner wedge
		local front = Vector3.new(0, halfSize.Z, halfSize.Y).Unit
		local right = Vector3.new(-halfSize.Y, halfSize.X, 0).Unit
		
		table.insert(axes, cframe.Rotation * front)
		table.insert(axes, cframe.Rotation * right)
	elseif shape == Enum.PartType.Cylinder then
		local horiz = (cframe.Position - part.Position) * Vector3.new(1, 0, 1)
		table.insert(axes, horiz.Unit)
	end
	
	return pA, utilGenericSat(
		utilCylinderSupport(cylinder),
		supportB,
		velocity, accel, dt,
		axes
	)
end

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

local function utilSolveActive(cylinder: Cylinder, planes: { Contact })
	if #planes == 0 then
		-- Nothing
		return
	elseif #planes == 1 then
		local a = planes[1]

		-- Push away from nA
		local nA = a.normal
		local d0, v0 = nA:Dot(a.mtv), nA:Dot(cylinder.velocity)

		cylinder.velocity -= nA * v0
		cylinder.position += nA * d0
	elseif #planes == 2 then
		local a, b = planes[1], planes[2]

		-- The only real direction we can move here is orthogonal to nA and nB
		local nA, nB = a.normal, b.normal
		local nAxB = nA:Cross(nB)

		local det = nAxB:Dot(nAxB)
		if math.abs(det) < 1e-6 then return end

		local d0, v0 = nA:Dot(a.mtv), nA:Dot(cylinder.velocity)
		local d1, v1 = nB:Dot(b.mtv), nB:Dot(cylinder.velocity)

		cylinder.position += (nB:Cross(nAxB) * d0 + nAxB:Cross(nA) * d1) / det
		cylinder.velocity -= (nB:Cross(nAxB) * v0 + nAxB:Cross(nA) * v1) / det
	else
		-- Concatenate all of the normal vectors
		-- as a column vector and solve by inverting the matrix 
		local a, b, c = planes[1], planes[2], planes[3]
		local nA, nB, nC = a.normal, b.normal, c.normal

		local det = nC:Dot(nA:Cross(nB))
		if math.abs(det) < 1e-6 then return end

		local d0, d1, d2 = nA:Dot(a.mtv), nB:Dot(b.mtv), nC:Dot(c.mtv)
		local v0, v1, v2 = nA:Dot(cylinder.velocity), nB:Dot(cylinder.velocity), nC:Dot(cylinder.velocity)
		cylinder.position += (d0 * nB:Cross(nC) + d1 * nC:Cross(nA) + d2 * nA:Cross(nB)) / det
		cylinder.velocity -= (v0 * nB:Cross(nC) + v1 * nC:Cross(nA) + v2 * nA:Cross(nB)) / det
	end
end


local function utilRelevantPlanes(cylinder: Cylinder, planes: { Contact }): { Contact }
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
		if plane.normal:Dot(cylinder.velocity) >= 0 then continue end
		if plane.normal:Dot(plane.mtv) <= 0.0 then continue end
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

function Cylinder.getContacts(cylinder: Cylinder, position: Vector3, size: Vector3): { BasePart }
	-- Stuff that would be considered as potential
	-- cylinders
	local primitives = {} :: { BasePart }

	for i, part in workspace:GetPartBoundsInBox(CFrame.new(position), size, cylinder.overlap) do
		if not part.CanCollide then
			continue
		elseif part:IsA("MeshPart") or part:IsA("UnionOperation") then
			continue
		else
			-- Otherwise, we can use the analytical solver
			table.insert(primitives, part :: BasePart)
		end
	end

	return primitives
end

function Cylinder.push(cylinder: Cylinder, what: BasePart, mtv: Vector3, linFac: number)
	if mtv:FuzzyEq(Vector3.zero) then return end

	-- We assume a typical PGS formulation.
	-- Also, this is exactly how you'd do a contact constraint
	-- in a typical physics engine :3
	local pb = what:GetClosestPointOnSurface(cylinder.position)
	local rA, rB = pb - cylinder.position, pb - what.Position

	-- Calculate the friction directions
	local normal = mtv.Unit
	local up = math.abs(normal.Y) < 0.99 and Vector3.yAxis or Vector3.xAxis
	local tangent = up:Cross(normal).Unit
	local bitangent = normal:Cross(tangent).Unit

	-- Get the current state variables
	local vB, wB = what.AssemblyLinearVelocity, what.AssemblyAngularVelocity
	local vA, wA = cylinder.velocity, cylinder.angularVelocity

	local cframeA, cframeB = CFrame.new(cylinder.position) * cylinder.rotation, what.CFrame
	local sA, sB = Vector3.new(cylinder.radius, cylinder.height, cylinder.radius) * 2, what.Size
	local iA, iB = utilInverseInertia(cylinder.rotation, sA, cylinder.mass), utilInverseInertia(cframeB.Rotation, sB, what.AssemblyMass)
	local mA, mB = 1.0 / cylinder.mass, 1.0 / what.AssemblyMass
	
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
	local maxFriction = math.max(-nDot * invEffMass, 0.0) * cylinder.friction

	-- Make a simplistic friction square (since, we don't need to do anything fancy
	-- for a rotation effect, lol.) Ideally, you'd use a 2x2 matrix inverse for this.
	local tImpulse = math.clamp(tDot / effMassT, -maxFriction, maxFriction)
	local bImpulse = math.clamp(bDot / effMassB, -maxFriction, maxFriction)
	local friction = tangent * tImpulse + bitangent * bImpulse

	-- Apply impulse on cylinder and candidate. Do not apply friction
	-- onto the cylinder's linear velocity (but angular velocity is okay.)
	local impulse = math.max(-nDot * invEffMass, 0.0) * normal + friction
	cylinder.angularVelocity += iA * rA:Cross(impulse)

	cylinder.position += mtv * linFac
	cylinder.velocity += mA * (impulse - friction) * linFac

	if not what.Anchored then
		what.AssemblyLinearVelocity -= impulse * mB
		what.AssemblyAngularVelocity -= iB * rB:Cross(impulse)
	end
end

function Cylinder.collide(c: Cylinder, dt: number)
	-- Get our contacts first
	local delta = c.velocity * dt + c.acceleration * dt * dt * 0.5
	local size = Vector3.new(c.radius, c.height, c.radius) * 2	
	local primitiveContacts = Cylinder.getContacts(
		c,
		c.position + delta,
		size + delta:Abs() * 2 + Vector3.one * 2
	)

	table.clear(c.contacts)	

	-- One CCD check is sufficient
	local ccd, ccdNormal, ccdPart, ccdPoint = dt, Vector3.zero, nil, Vector3.zero
	for i, part in primitiveContacts do
		local contact, _, _, toiNormal, toi = utilCylinderCastPrimitive(c, part, c.velocity, c.acceleration, dt)
		if toi > 0.0 and toi < ccd then
			ccdNormal = toiNormal
			ccdPoint = contact
			ccdPart = part
			ccd = toi
		end	
	end

	c.position += c.velocity * ccd + 0.5 * c.acceleration * ccd * ccd
	if ccd < dt then
		-- Minimize dot(disp, normal), given velocity + normal*lambda_1
		local dot = math.min(c.velocity:Dot(ccdNormal), 0.0)
		c.velocity -= ccdNormal * dot
		c.position += c.velocity * (dt - ccd) + 0.5 * c.acceleration * (dt - ccd)^2
		dt -= ccd

		table.insert(c.contacts, {
			position = ccdPoint,
			normal = ccdNormal,
			mtv = Vector3.zero,
			part = ccdPart :: BasePart,
		})
	end
	
	for i, part in primitiveContacts do
		local contact, normal, push, _, toi = utilCylinderCastPrimitive(c, part, c.velocity, c.acceleration, dt)
		if toi > 0.0 then continue end

		-- Clip velocity and use the current MTV
		Cylinder.push(c, part, normal * push, part.Anchored and 0.0 or 1.0)
		table.insert(c.contacts, {
			position = contact,
			normal = normal,
			mtv = push * normal,
			part = part
		})
	end

	utilSolveActive(c, utilRelevantPlanes(c, c.contacts))	
end

Cylinder.utilInverseInertia = utilInverseInertia

return Cylinder