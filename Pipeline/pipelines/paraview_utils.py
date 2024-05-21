from paraview.simple import *
sphereInstance = Sphere()
sphereInstance.Radius = 1.0
sphereInstance.Center[1] = 1.0
sphereDisplay = Show(sphereInstance)
view = Render()
sphereDisplay.Opacity = 0.5
Render(view)