#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
a144_edgeTrue_subnet_and_genevtu = XMLUnstructuredGridReader(FileName=['/data/SonData/BlenderExercises/vtk_test/data/layout/paraview_1.4.4/1.4.4_edgeTrue_subnet_and_gene.vtu'])
a144_edgeTrue_subnet_and_genevtu.PointArrayStatus = ['Degree', 'SUID', 'X1', 'X2', 'change_direction', 'degree', 'diff', 'enrichment_node_smooth_score', 'enrichment_subnet_score', 'enrichment_subnet_score_with_sign', 'level', 'n_down', 'n_genes', 'n_up', 'nbh_total_score', 'net_type', 'node_smooth_score', 'object_type', 'pfdr', 'pval', 'radius', 'score', 'selected', 'shared.name', 'subnet_pvalue', 'subnet_qvalue', 'subnet_score', 'weight']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1399, 838]

# show data in view
a144_edgeTrue_subnet_and_genevtuDisplay = Show(a144_edgeTrue_subnet_and_genevtu, renderView1)


# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(a144_edgeTrue_subnet_and_genevtuDisplay, ('POINTS', 'enrichment_subnet_score'))

# rescale color and/or opacity maps used to include current data range
a144_edgeTrue_subnet_and_genevtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
a144_edgeTrue_subnet_and_genevtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'enrichment_subnet_score'
enrichment_subnet_scoreLUT = GetColorTransferFunction('enrichment_subnet_score')

# Properties modified on renderView1
renderView1.Background = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [0.04854755476117134, 0.25337615609169006, 0.028104693815075764]
renderView1.CameraFocalPoint = [0.04854755476117134, 0.25337615609169006, 0.0]
renderView1.CameraParallelScale = 0.00727403001611663

# save screenshot
SaveScreenshot('/data/SonData/BlenderExercises/vtk_test/data/layout/paraview_1.4.4/subnet_network.png', renderView1, ImageResolution=[13990, 8380],
    ImageQuality=80)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.04812671422904469, 0.25341441432188333, 0.016798030226780922]
renderView1.CameraFocalPoint = [0.04812671422904469, 0.25341441432188333, 0.0]
renderView1.CameraParallelScale = 0.00727403001611663

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).