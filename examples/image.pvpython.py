#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
imageicefloes = XMLUnstructuredGridReader(FileName=['./image/image.icefloes.1.vtu', './image/image.icefloes.2.vtu', './image/image.icefloes.3.vtu', './image/image.icefloes.4.vtu', './image/image.icefloes.5.vtu', './image/image.icefloes.6.vtu', './image/image.icefloes.7.vtu', './image/image.icefloes.8.vtu', './image/image.icefloes.9.vtu', './image/image.icefloes.10.vtu', './image/image.icefloes.11.vtu', './image/image.icefloes.12.vtu', './image/image.icefloes.13.vtu', './image/image.icefloes.14.vtu', './image/image.icefloes.15.vtu', './image/image.icefloes.16.vtu', './image/image.icefloes.17.vtu', './image/image.icefloes.18.vtu', './image/image.icefloes.19.vtu', './image/image.icefloes.20.vtu', './image/image.icefloes.21.vtu', './image/image.icefloes.22.vtu', './image/image.icefloes.23.vtu', './image/image.icefloes.24.vtu', './image/image.icefloes.25.vtu', './image/image.icefloes.26.vtu', './image/image.icefloes.27.vtu', './image/image.icefloes.28.vtu', './image/image.icefloes.29.vtu', './image/image.icefloes.30.vtu', './image/image.icefloes.31.vtu', './image/image.icefloes.32.vtu', './image/image.icefloes.33.vtu', './image/image.icefloes.34.vtu', './image/image.icefloes.35.vtu', './image/image.icefloes.36.vtu', './image/image.icefloes.37.vtu', './image/image.icefloes.38.vtu', './image/image.icefloes.39.vtu', './image/image.icefloes.40.vtu', './image/image.icefloes.41.vtu', './image/image.icefloes.42.vtu', './image/image.icefloes.43.vtu', './image/image.icefloes.44.vtu', './image/image.icefloes.45.vtu', './image/image.icefloes.46.vtu', './image/image.icefloes.47.vtu', './image/image.icefloes.48.vtu', './image/image.icefloes.49.vtu', './image/image.icefloes.50.vtu'])
imageicefloes.PointArrayStatus = ['Density [kg m^-3]', 'Thickness [m]', 'Diameter (contact) [m]', 'Diameter (areal) [m]', 'Circumreference  [m]', 'Horizontal surface area [m^2]', 'Side surface area [m^2]', 'Volume [m^3]', 'Mass [kg]', 'Moment of inertia [kg m^2]', 'Linear velocity [m s^-1]', 'Linear acceleration [m s^-2]', 'Sum of forces [N]', 'Angular position [rad]', 'Angular velocity [rad s^-1]', 'Angular acceleration [rad s^-2]', 'Sum of torques [N*m]', 'Fixed in space [-]', 'Free to rotate [-]', 'Enabled [-]', 'Contact stiffness (normal) [N m^-1]', 'Contact stiffness (tangential) [N m^-1]', 'Contact viscosity (normal) [N m^-1 s]', 'Contact viscosity (tangential) [N m^-1 s]', 'Contact friction (static) [-]', 'Contact friction (dynamic) [-]', "Young's modulus [Pa]", "Poisson's ratio [-]", 'Tensile strength [Pa]', 'Compressive strength prefactor [m^0.5 Pa]', 'Ocean drag coefficient (vertical) [-]', 'Ocean drag coefficient (horizontal) [-]', 'Atmosphere drag coefficient (vertical) [-]', 'Atmosphere drag coefficient (horizontal) [-]', 'Contact pressure [Pa]', 'Number of contacts [-]', 'Granular stress [Pa]', 'Ocean stress [Pa]', 'Atmosphere stress [Pa]']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2478, 1570]

# show data in view
imageicefloesDisplay = Show(imageicefloes, renderView1)
# trace defaults for the display properties.
imageicefloesDisplay.Representation = 'Surface'
imageicefloesDisplay.AmbientColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.ColorArrayName = [None, '']
imageicefloesDisplay.OSPRayScaleArray = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
imageicefloesDisplay.SelectOrientationVectors = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.ScaleFactor = 6.050000000000001
imageicefloesDisplay.SelectScaleArray = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.GlyphType = 'Arrow'
imageicefloesDisplay.GlyphTableIndexArray = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.DataAxesGrid = 'GridAxesRepresentation'
imageicefloesDisplay.PolarAxes = 'PolarAxesRepresentation'
imageicefloesDisplay.ScalarOpacityUnitDistance = 64.20669746996803
imageicefloesDisplay.GaussianRadius = 3.0250000000000004
imageicefloesDisplay.SetScaleArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
imageicefloesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
imageicefloesDisplay.OpacityArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
imageicefloesDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
imageicefloesDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
imageicefloesDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [30.5, 11.0, 10000.0]
#renderView1.CameraFocalPoint = [30.5, 11.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(Input=imageicefloes,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
glyph1.Vectors = ['POINTS', 'Angular acceleration [rad s^-2]']
glyph1.ScaleFactor = 6.050000000000001
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.Scalars = ['POINTS', 'Diameter (areal) [m]']
glyph1.Vectors = ['POINTS', 'Angular position [rad]']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 1.0
glyph1.GlyphMode = 'All Points'

# get color transfer function/color map for 'Diameterarealm'
diameterarealmLUT = GetColorTransferFunction('Diameterarealm')

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.LookupTable = diameterarealmLUT
glyph1Display.OSPRayScaleArray = 'Diameter (areal) [m]'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 6.1000000000000005
glyph1Display.SelectScaleArray = 'Diameter (areal) [m]'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Diameter (areal) [m]'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 3.0500000000000003
glyph1Display.SetScaleArray = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'

# update the view to ensure updated data information
renderView1.Update()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# rescale color and/or opacity maps used to exactly fit the current data range
glyph1Display.RescaleTransferFunctionToDataRange(False, True)

# Rescale transfer function
diameterarealmLUT.RescaleTransferFunction(0.0, 5.0)

# get opacity transfer function/opacity map for 'Diameterarealm'
diameterarealmPWF = GetOpacityTransferFunction('Diameterarealm')

# Rescale transfer function
diameterarealmPWF.RescaleTransferFunction(0.0, 5.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
diameterarealmLUT.ApplyPreset('Black, Blue and White', True)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [30.75, 10.99999962002039, 125.01319337485243]
#renderView1.CameraFocalPoint = [30.75, 10.99999962002039, 0.0]
#renderView1.CameraParallelScale = 20.56227912200389

# save animation
#SaveAnimation('./image/image.avi', renderView1, ImageResolution=[1239, 785],
#    FrameRate=10,
#    FrameWindow=[0, 49])

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [30.75, 10.99999962002039, 125.01319337485243]
#renderView1.CameraFocalPoint = [30.75, 10.99999962002039, 0.0]
#renderView1.CameraParallelScale = 20.56227912200389

# save animation
SaveAnimation('./image/image.png', renderView1, ImageResolution=[1239, 785],
    FrameRate=10,
    FrameWindow=[0, 49])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [30.75, 10.99999962002039, 125.01319337485243]
#renderView1.CameraFocalPoint = [30.75, 10.99999962002039, 0.0]
#renderView1.CameraParallelScale = 20.56227912200389

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
