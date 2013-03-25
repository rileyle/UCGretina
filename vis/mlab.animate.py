from mayavi import mlab

@mlab.animate(delay=100)

def anim():
    f = mlab.gcf()
    while(1):
        # Aim the camera at the center of the target 
        # without resetting "camera up"
        mlab.view(focalpoint=(0,0,0),distance=4000.0,reset_roll=False)

        # Rotate the camera about "camera up" (angle in degrees)
        f.scene.camera.azimuth(5)

        f.scene.render()
        yield

a = anim() # Starts the animation.
