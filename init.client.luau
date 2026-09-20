--!strict

-- Quick bootstrapper to test the movement system backend
-- magicoal_nerb

local Quakem = require("@game/ReplicatedStorage/quakem/Quakem")

local RunService = game:GetService("RunService")
local Players = game:GetService("Players")

local LocalPlayer = assert(Players.LocalPlayer)
local CurrentCamera = assert(workspace.CurrentCamera)

local function characterAdded()
	local character = assert(LocalPlayer.Character)
	local humanoid = character:WaitForChild("Humanoid") :: Humanoid
	
	local quake = Quakem.new(humanoid)
	humanoid.EvaluateStateMachine = false
	
	RunService:BindToRenderStep("quakem-update", Enum.RenderPriority.Camera.Value - 1, function(dt: number)
		dt = math.min(dt, 1/30)
		if UserSettings().GameSettings.RotationType == Enum.RotationType.CameraRelative then
			local lookVector = CurrentCamera.CFrame.ZVector
			quake.rotation = CFrame.Angles(0, math.atan2(lookVector.X, lookVector.Z), 0)
			quake.angularVelocity = Vector3.zero
		else
			Quakem.think(quake, dt)
		end
		
		Quakem.step(quake, dt)
	end)
end

if LocalPlayer.Character then
	characterAdded()
end

LocalPlayer.CharacterAdded:Connect(characterAdded)