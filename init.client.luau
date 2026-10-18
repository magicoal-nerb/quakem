local Quakem = require(script.Parent:WaitForChild("quakem"):WaitForChild("Quakem"))

-- Example code!! This is not deterministic on its own!!
-- magicoal_nerb :3

local StarterGui = game:GetService("StarterGui")
local RunService = game:GetService("RunService")
local Players = game:GetService("Players")

local LocalPlayer = assert(Players.LocalPlayer)
local CurrentCamera = assert(workspace.CurrentCamera)

local ResetBindable = Instance.new("BindableEvent")
local Event: RBXScriptConnection

local function characterAdded()
	local character = assert(LocalPlayer.Character)
	local humanoid = character:WaitForChild("Humanoid") :: Humanoid
	local root = character:WaitForChild("HumanoidRootPart")
	
	local quake = Quakem.new(humanoid)
	humanoid.EvaluateStateMachine = false
	humanoid.JumpPower = 25
	
	RunService:BindToRenderStep("quakem-update", Enum.RenderPriority.Camera.Value - 1, function(dt: number)
		if UserSettings().GameSettings.RotationType == Enum.RotationType.CameraRelative then
			local lookVector = CurrentCamera.CFrame.ZVector
			quake.rotation = CFrame.Angles(0, math.atan2(lookVector.X, lookVector.Z), 0)
			quake.angularVelocity = Vector3.zero
		else
			Quakem.think(quake, dt)
		end
		
		Quakem.step(quake, dt)
	end)
	
	Event = ResetBindable.Event:Connect(function()
		quake.position = workspace.SpawnLocation.Position
		quake.velocity = Vector3.zero
	end)
end

if LocalPlayer.Character then
	characterAdded()
end

LocalPlayer.CharacterAdded:Connect(characterAdded)

while task.wait(0.5) do
	local succ, err = pcall(StarterGui.SetCore, StarterGui, "ResetButtonCallback", ResetBindable)
	if succ then
		break
	end
end