#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
logoicefloes = XMLUnstructuredGridReader(FileName=['./logo/logo.icefloes.1.vtu', './logo/logo.icefloes.2.vtu', './logo/logo.icefloes.3.vtu', './logo/logo.icefloes.4.vtu', './logo/logo.icefloes.5.vtu', './logo/logo.icefloes.6.vtu', './logo/logo.icefloes.7.vtu', './logo/logo.icefloes.8.vtu', './logo/logo.icefloes.9.vtu', './logo/logo.icefloes.10.vtu', './logo/logo.icefloes.11.vtu', './logo/logo.icefloes.12.vtu', './logo/logo.icefloes.13.vtu', './logo/logo.icefloes.14.vtu', './logo/logo.icefloes.15.vtu', './logo/logo.icefloes.16.vtu', './logo/logo.icefloes.17.vtu', './logo/logo.icefloes.18.vtu', './logo/logo.icefloes.19.vtu', './logo/logo.icefloes.20.vtu', './logo/logo.icefloes.21.vtu', './logo/logo.icefloes.22.vtu', './logo/logo.icefloes.23.vtu', './logo/logo.icefloes.24.vtu', './logo/logo.icefloes.25.vtu', './logo/logo.icefloes.26.vtu', './logo/logo.icefloes.27.vtu', './logo/logo.icefloes.28.vtu', './logo/logo.icefloes.29.vtu', './logo/logo.icefloes.30.vtu', './logo/logo.icefloes.31.vtu', './logo/logo.icefloes.32.vtu', './logo/logo.icefloes.33.vtu', './logo/logo.icefloes.34.vtu', './logo/logo.icefloes.35.vtu', './logo/logo.icefloes.36.vtu', './logo/logo.icefloes.37.vtu', './logo/logo.icefloes.38.vtu', './logo/logo.icefloes.39.vtu', './logo/logo.icefloes.40.vtu', './logo/logo.icefloes.41.vtu', './logo/logo.icefloes.42.vtu', './logo/logo.icefloes.43.vtu', './logo/logo.icefloes.44.vtu', './logo/logo.icefloes.45.vtu', './logo/logo.icefloes.46.vtu', './logo/logo.icefloes.47.vtu', './logo/logo.icefloes.48.vtu', './logo/logo.icefloes.49.vtu', './logo/logo.icefloes.50.vtu'])
logoicefloes.PointArrayStatus = ['Density [kg m^-3]', 'Thickness [m]', 'Diameter (contact) [m]', 'Diameter (areal) [m]', 'Circumreference  [m]', 'Horizontal surface area [m^2]', 'Side surface area [m^2]', 'Volume [m^3]', 'Mass [kg]', 'Moment of inertia [kg m^2]', 'Linear velocity [m s^-1]', 'Linear acceleration [m s^-2]', 'Sum of forces [N]', 'Angular position [rad]', 'Angular velocity [rad s^-1]', 'Angular acceleration [rad s^-2]', 'Sum of torques [N*m]', 'Fixed in space [-]', 'Free to rotate [-]', 'Enabled [-]', 'Contact stiffness (normal) [N m^-1]', 'Contact stiffness (tangential) [N m^-1]', 'Contact viscosity (normal) [N m^-1 s]', 'Contact viscosity (tangential) [N m^-1 s]', 'Contact friction (static) [-]', 'Contact friction (dynamic) [-]', "Young's modulus [Pa]", "Poisson's ratio [-]", 'Tensile strength [Pa]', 'Compressive strength prefactor [m^0.5 Pa]', 'Ocean drag coefficient (vertical) [-]', 'Ocean drag coefficient (horizontal) [-]', 'Atmosphere drag coefficient (vertical) [-]', 'Atmosphere drag coefficient (horizontal) [-]', 'Contact pressure [Pa]', 'Number of contacts [-]', 'Granular stress [Pa]', 'Ocean stress [Pa]', 'Atmosphere stress [Pa]']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2478, 1570]

# show data in view
logoicefloesDisplay = Show(logoicefloes, renderView1)
# trace defaults for the display properties.
logoicefloesDisplay.Representation = 'Surface'
logoicefloesDisplay.AmbientColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.ColorArrayName = [None, '']
logoicefloesDisplay.OSPRayScaleArray = 'Angular acceleration [rad s^-2]'
logoicefloesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
logoicefloesDisplay.SelectOrientationVectors = 'Angular acceleration [rad s^-2]'
logoicefloesDisplay.ScaleFactor = 6.050000000000001
logoicefloesDisplay.SelectScaleArray = 'Angular acceleration [rad s^-2]'
logoicefloesDisplay.GlyphType = 'Arrow'
logoicefloesDisplay.GlyphTableIndexArray = 'Angular acceleration [rad s^-2]'
logoicefloesDisplay.DataAxesGrid = 'GridAxesRepresentation'
logoicefloesDisplay.PolarAxes = 'PolarAxesRepresentation'
logoicefloesDisplay.ScalarOpacityUnitDistance = 64.20669746996803
logoicefloesDisplay.GaussianRadius = 3.0250000000000004
logoicefloesDisplay.SetScaleArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
logoicefloesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
logoicefloesDisplay.OpacityArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
logoicefloesDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
logoicefloesDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
logoicefloesDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
logoicefloesDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [30.5, 11.0, 10000.0]
#renderView1.CameraFocalPoint = [30.5, 11.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(Input=logoicefloes,
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
#SaveAnimation('./logo/logo.avi', renderView1, ImageResolution=[1239, 785],
#    FrameRate=10,
#    FrameWindow=[0, 49])

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [30.75, 10.99999962002039, 125.01319337485243]
#renderView1.CameraFocalPoint = [30.75, 10.99999962002039, 0.0]
#renderView1.CameraParallelScale = 20.56227912200389

# save animation
SaveAnimation('./logo/logo.png', renderView1, ImageResolution=[1239, 785],
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
