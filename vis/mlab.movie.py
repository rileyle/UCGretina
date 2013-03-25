from mayavi import mlab

@mlab.animate(delay=100)           # (delay between frames in ms)

def anim():
    f = mlab.gcf()
    frameNumber = 0
    for i in range(0,360,5):
        # Aim the camera at the center of the target 
        # without resetting "camera up"
        mlab.view(focalpoint=(0,0,0),distance=4000.0,reset_roll=False)

        # Rotate the camera about "camera up" (angle in degrees)
        f.scene.camera.azimuth(5)

        f.scene.render()
        yield

#        filename = 'frame{}'.format(i)
        frameNumber += 1
        filename = 'frame{:03d}'.format(frameNumber)
        filename += '.png'
        print filename
        mlab.savefig(filename)

a = anim() # Starts the animation.
