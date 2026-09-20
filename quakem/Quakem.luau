--!strict

-- main quake movement system :P
-- to change between colliders, set the Collider variable to capsule/cylinder
-- magicoal_nerb

local Collider = require("./Capsule")
local Quakem = {}

const kVecXZ = Vector3.new(1, 0, 1)
const kFloorNormal = math.sqrt2 / 2
const kGroundAcceleration = 220
const kGroundFriction = 14
const kAirAcceleration = 600
const kMaxAirSpeed = 2.4
const kGravity = 90

const kTurnStiffness = 50.0
const kTurnSpeed = 8.0

-- Visual effect
const kWallFriction = 0.5

export type Quakem = {
	position: Vector3,
	angularVelocity: Vector3,
	acceleration: Vector3,
	velocity: Vector3,

	humanoid: Humanoid,
	rotation: CFrame,
	friction: number,
	mass: number,

	height: number,
	radius: number,
	overlap: OverlapParams,
	params: RaycastParams,
	rootPart: BasePart,

	contacts: { Collider.Contact },
}

local function utilPlayerHeight(humanoid: Humanoid): number
	if humanoid.RigType == Enum.HumanoidRigType.R6 then
		-- If it is R6, we just use the default size
		return 3
	else
		-- Otherwise, we will just use the height provided
		-- by the humanoid
		return humanoid.HipHeight + 0.9
	end
end

local function utilAngleBetween(a: Vector3, b: Vector3): number
	-- sin = cross(a, b).y, cos = dot(a, b)
	local sin = a.X * b.Z - a.Z * b.X
	local cos = a.X * b.X + a.Z * b.Z

	return math.atan2(sin, cos)
end

local function utilAccelerate(
	input: Vector3,
	velocity: Vector3,
	accelerate: number,
	maxVelocity: number,
	dt: number
): Vector3
	-- Constrain min(dot(v, i + i*dt), maxVel)
	local projection: number = velocity:Dot(input)
	local accel: number = math.max(math.min(accelerate * dt, maxVelocity - projection), 0.0)
	return velocity + accel * input
end

function Quakem.new(humanoid: Humanoid): Quakem
	local rootPart = assert(humanoid.RootPart, "quakem(): no rootpart!")
	local overlap = OverlapParams.new()
	local params = RaycastParams.new()
	params.FilterDescendantsInstances = { rootPart.Parent :: BasePart }
	params.FilterType = Enum.RaycastFilterType.Exclude
	params.RespectCanCollide = true
	
	overlap.FilterDescendantsInstances = { rootPart.Parent :: BasePart }
	overlap.FilterType = Enum.RaycastFilterType.Exclude
	overlap.RespectCanCollide = true
	
	return {
		rootPart = rootPart,

		position = rootPart.Position,
		angularVelocity = rootPart.AssemblyAngularVelocity,
		velocity = rootPart.AssemblyLinearVelocity,

		humanoid = humanoid,
		rotation = CFrame.identity,
		mass = rootPart.AssemblyMass,

		height = utilPlayerHeight(humanoid),
		radius = 1,
		friction = kWallFriction,
		
		overlap = overlap,
		params = params,
		
		contacts = {}
	} :: Quakem
end

function Quakem.think(q: Quakem, dt: number)
	local wishDir, angularVelocity = q.humanoid.MoveDirection, q.angularVelocity
	local desiredAngular = if wishDir:Dot(wishDir) > 0.0
		then Vector3.yAxis * kTurnSpeed * utilAngleBetween(-wishDir, q.rotation.ZVector)
		else Vector3.zero

	local floor, floorPart = Vector3.zero, nil
	for i, contact in q.contacts do
		local normal = contact.normal
		if normal.Y > kFloorNormal and floor.Y < normal.Y then
			floorPart = contact.part
			floor = normal
		end
	end

	if floorPart then
		-- If on the floor, try to align yourself to the angular
		-- velocity of the ground
		desiredAngular += floorPart.AssemblyAngularVelocity
	else
		local projection = q.velocity * kVecXZ
		if projection:Dot(projection) > 0.1 then
			-- Look towards our desired velocity when in midair
			desiredAngular = Vector3.yAxis * kTurnSpeed * utilAngleBetween(-projection.Unit, q.rotation.ZVector)
		end
	end

	-- Implicit: inertia * Vector3.Up * kP * (desiredAngular + floorAngular - angularVelocity.Y)
	-- Let w_t = w_d + w_f
	-- I dw/dt = I*j*(kP * (w_t - w(0)))
	-- dw_y/dt = kP * (w_t - w(0))
	-- dw_dt/dt + kP*w(0) = kP*w_t
	-- u(t) = e^int kP dt <=> e^(kP*t)
	-- <=> e^(kP*dt)w_y = e^(kP*t)*kP*w_t
	-- <=> e^(kP*dt)w_y(t) - w_y = w_t(e^(kP*t) - 1)
	-- <=> e^(-kP*dt) w_y + (1-e^(-kP*dt))*w_t
	local kPt = math.exp(-kTurnStiffness * dt)
	local invCapsuleInertia = Collider.utilInverseInertia(
		q.rotation, Vector3.new(q.radius, q.height, q.radius) * 2.0,
		q.mass
	)

	q.angularVelocity += (kPt - 1.0) * q.angularVelocity + (1.0 - kPt) * desiredAngular * Vector3.yAxis
	q.angularVelocity *= Vector3.yAxis
	q.rotation = CFrame.Angles(0, q.angularVelocity.Y * dt, 0) * q.rotation
end

function Quakem.step(q: Quakem, dt: number)
	local wishDirection = q.humanoid.MoveDirection

	local floor, floorPart = Vector3.zero, nil
	for i, contact in q.contacts do
		local normal = contact.normal
		if normal.Y > kFloorNormal and floor.Y < normal.Y then
			floorPart = contact.part
			floor = normal
		end
	end
	
	debug.profilebegin("quakem::collide")

	local state: Enum.HumanoidStateType
	if q.humanoid.Jump and floor.Y ~= 0.0 then
		-- Just jump!
		state = Enum.HumanoidStateType.Jumping
		q.acceleration = Vector3.zero
		q.velocity = Vector3.new(q.velocity.X, 25, q.velocity.Z)
	elseif floor.Y ~= 0.0 and floorPart then
		-- Then, we're on the ground. Just apply basic
		-- friction and go.
		local v, w = floorPart.AssemblyLinearVelocity, floorPart.AssemblyAngularVelocity
		local floorVelocity = v + w:Cross(q.position - floorPart.Position)

		-- Apply a friction force (that goes towards floorVelocity)
		-- dv/dt = k(v_t - v) => dv/(v_t - v) = kdt
		-- -ln(v_t - v) = kt + C
		-- v_t - v = Ce^-kt
		-- Know: v(0) = C <=> v_t - v
		-- v = v_t - (v_t - v)e^-kt
		local exponential = math.exp(-kGroundFriction * dt)
		q.velocity = floorVelocity - (floorVelocity - q.velocity) * exponential
		q.velocity = utilAccelerate(
			wishDirection,
			q.velocity * kVecXZ,
			kGroundAcceleration,
			q.humanoid.WalkSpeed,
			dt
		) * kVecXZ

		state = Enum.HumanoidStateType.Running
		q.acceleration = -floor:Dot(wishDirection * wishDirection:Dot(q.velocity)) * Vector3.yAxis
		q.velocity -= floor:Dot(q.velocity) * Vector3.yAxis + floor
	else
		-- Then we're airborne
		local projection = utilAccelerate(
			wishDirection,
			q.velocity * kVecXZ,
			kAirAcceleration,
			kMaxAirSpeed,
			dt
		)

		state = Enum.HumanoidStateType.Freefall
		q.acceleration = -kGravity * Vector3.yAxis
		q.velocity = (q.velocity * Vector3.yAxis) + projection
		q.velocity -= kGravity * Vector3.yAxis * dt
	end

	Collider.collide(q, dt)
	debug.profileend()
	
	q.humanoid:ChangeState(state)
	q.rootPart.CFrame = CFrame.new(q.position) * q.rotation
	q.rootPart.AssemblyLinearVelocity = q.velocity
	q.rootPart.AssemblyAngularVelocity = q.angularVelocity
end

return Quakem