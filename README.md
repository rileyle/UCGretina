# UCGretina #

Cite: [L.A.Riley, D.Weisshaar, H.L.Crawford et al., UCGretina GEANT4 simulation of the GRETINA Gamma-Ray Energy Tracking Array, Nucl. Instr. Meth. A1003, 165305 (2021)](https://doi.org/10.1016/j.nima.2021.165305)

## Compile and install ##

Install version [Geant4-10.7.4 of the Geant4 libraries](https://geant4.web.cern.ch/download/all). You will need the data files for low energy electromagnetic processes, photon evaporation, and radioactive decay.

The model of the GRETINA scanning table uses version 2.0.3 of the external [CADMesh](https://github.com/christopherpoole/cadmesh) package. 

Set up your environment (consider adding this to your `.bashrc`):

    $ source <Path to Geant4>/bin/geant4.sh
    $ source <Path to Geant4>/share/Geant4-10.7.4/geant4make/geant4make.sh

Compile:

    $ make

To use the LBL scanning table:

    $ make SCANNING=1

(produces the binary UCGretina_Scan)

To use the liquid hydrogen target:

    $ make LHTARGET=1

(produces the binary UCGretina_LH)

To include nuclear polarization (alignment) of the reaction product in the `Reaction` class:

    $ make POL=1

(produces the binary `UCGretina_Pol`)
Implementation and validation of this capability is described here: [C. Morse, H. L. Crawford, A. O. Macchiavelli et al., The polarization sensitivity of GRETINA, Nucl. Instr. Meth. A1025, 166155 (2022)](https://doi.org/10.1016/j.nima.2021.166155)

To activate neutron-related processes in the physics list (required for the `neutron` source type:

    $ make NEUTRONS=1

(does not affect the executable name, can be used with any of the above flags)

Executables are automatically installed in

    $G4WORKDIR/bin/$G4SYSTEM

(which is added to your path when you source `geant4make.sh`)

## Examples ##

Several examples are in the examples subdirectory, including
illustrations of fitting simulations to measured source and in-beam
spectra. Makefiles are provided in the examples for sorting simulated
mode 2 output with the
[GRUTinizer](https://github.com/pcbend/GRUTinizer) and
[GrROOT](https://github.com/wimmer-k/GrROOT)
packages.

## Selected Macro File Commands ##

### Geometry ###

Five text files present in the working directory in which `UCGretina` is run determine the geometry of the crystals (`asolid`), including coaxial and back passive layers, the segmentation of crystals for readout (`aslice`), the assembly of crystals into quads (`acluster`), the aluminum vacuum jacket surrounding the crystals in each quad (`awalls`), and the placement of modules into the array (`aeuler`). The detector construction code is a modified version of the `AgataDetectorArray` class in the AGATA simulation code. The file `./GretinaGeometry/geometry_description_agata` describes the geometry file formats.

Geometry files for several array configurations are provided in the `./GretinaGeometry` directory. The bash script `change_geometry.sh` is provided to create soft links in the current working directory to the geometry files with the path and root name supplied as a command line argument. For example,

    $ ./change_geometry.sh <PATH TO GretinaGeometry>/GretinaNSCL/G120C4

creates links to the 7-quad configuration used in the GRETINA commissioning at the NSCL.

Optional commands for setting target parameters:

    /Target/Material <material name>
    
    /Target/X_length <double> <unit>
    
    /Target/Y_length <double> <unit>
    
    /Target/SetPosition_X <double> <unit>
    
    /Target/SetPosition_Y <double> <unit>
    
    /Target/SetPosition_Z <double> <unit>
    
    /Target/Thickness <double> <unit>
    
    /Target/ScaleDensity <double>
    
    /Target/Sled

Mandatory command for building a target/source frame:

    /Target/Construct

Optional commands for setting beam-tube geometry:

    /BeamTube/R_min <double> <unit>
    
    /BeamTube/R_max <double> <unit>
    
    /BeamTube/Length <double> <unit>

Mandatory command for building the beam tube:

    /BeamTube/Construct

Mandatory command for building the WU chamber:

    /WUChamber/Construct

Mandatory command for building the GRETA chamber (spherical chamber with 36 cm OD and 1.5 mm wall thickness, and a beam pipe with 4 cm OD and 1.5 mm wall thickness):

    /GretaChamber/Construct

Optional commands for setting GRETA chamber geometry:

    /GretaChamber/R_min <double> <unit>
    
    /GretaChamber/R_max <double> <unit>

Optional commands for including GRETINA-related passive material:

    /Gretina/detector/enableCapsules
    
    /Gretina/detector/enableCryostats
    
    /Gretina/Shell < full || north || south || Greta || GretaLH || Greta_North || Greta_South || GretaLH_North || GretaLH_South >

Optional commands setting the offset of the north and south halves of the mounting shell (positive values correspond to backing away from the target):

    /Gretina/NorthOffset <double> <unit>
    
    /Gretina/SouthOffset <double> <unit>

Optional command to omit the GRETINA detectors:

    /Gretina/NoDetectors

Optional command to include a model of the S800 quadrupole and gate valve:

    /Gretina/S800

### Gamma-Ray Angular Correlations (see also /examples/sources/co60) ###

Gamma-ray angular correlations are built into the `G4PhotonEvaporation/G4GammaTransition` classes (starting with geant4.10.4). This functionality is disabled by default but can be enabled with: 

    /PhysicsList/AngularCorrelations true
    /PhysicsList/SetGammaPolarization true

The `./examples/sources/co60` example includes a simulation of the angular correlations in the 4 -> 2 -> 0 cascade in <SUP>60</SUP>Ni and also shows how to simulate isotropic distributions for comparison with correlated ones.

### Initialize the Run Manager ###

After issuing the above commands as needed to pre-configure the simulation, the run manager must be initialized:

    /run/initialize

### In-beam Simulations ###

Mandatory command after all `/BeamOut/` and `/BeamIn/` commands:

    /BeamOut/Update

Mandatory commands related to the incoming beam:

    /BeamIn/A <A> /BeamIn/Z <Z>

> Mass number and atomic number

    /BeamIn/KEu <double> <unit>

> Kinetic energy per nucleon of the incoming beam

Optional commands related to the incoming beam:

    /BeamIn/Dpp <double>

> Momentum acceptance (dp/p) for the incoming beam. (Set this to 0 if you are providing the dta spectrum of the incoming beam using the `/BeamIn/dtaFile` command.)

    /BeamIn/dtaFile <filename>

> file name for the dta spectrum of the incoming beam. This is a text file with format: 

>       <Minimum DTA [%]> <Maximum DTA [%]> <DTA bin width [%]> 
>       <Bin 1 counts> 
>       <Bin 2 counts>
>       <Bin 3 counts>
>       ...

> If this command is present, the incoming beam energy is set by drawing randomly from the dta distribution centered on the beam energy specified by the `/BeamIn/KEu` command. (The dta spectrum specifies the momentum acceptance of the incoming beam, so the momentum acceptance parameter should be set to zero: `/BeamIn/Dpp 0`.)

    /BeamIn/Focus/X <double> <unit>
    
    /BeamIn/Focus/Y <double> <unit>

> Offsets of the emission point of the incoming beam. (Z defaults to -50 cm. If you change this, make sure it is upstream of the target!)

    /BeamIn/Focus/DX <double> <unit>
    
    /BeamIn/Focus/DY <double> <unit>

> Horizontal and vertical widths of the beam spot at the emission point (not on target)

    /BeamIn/Focus/Ata0 <double> <unit>
    
    /BeamIn/Focus/Bta0 <double> <unit>

> Direction of the incoming beam (dispersive and nondispersive angles, respectively)

    /BeamIn/Focus/maxAta <double> <unit> 
    
    /BeamIn/Focus/maxBta <double> <unit>

> Angular divergences of the incoming beam in the dispersive and nondispersive directions, respectively.

    /BeamIn/Excitation <double> <unit>

> Excitation energy of the incoming beam. This generally only makes sense in simulations of the the gamma decay of stationary excited nuclei, for which the `/BeamOut/Source` command is present. 

Mandatory commands related to the outgoing beam:

    /BeamOut/DA <int> 
    
    /BeamOut/DZ <int>

> Changes in mass number and atomic number of the reaction. The incoming beam has (A,Z) and the outgoing reaction product has (A+DA, Z+DZ)

    /BeamOut/TargetA <int>
    
    /BeamOut/TargetZ <int>

> Mass number and charge number of the target nucleus.

    /BeamOut/ProjectileExcitation <double> <unit>

> Excitation energy of the beam-like reaction product. (The target-like reaction product is not excited.) This energy must correspond to a state described in the level data file described below. 

    /BeamOut/TargetExcitation <double> <unit>

> Excitation energy of the target-like reaction product. (The beam-like reaction product is not excited.) This energy must correspond to a state described in the level data file described below. 

    /BeamOut/LevelDataFile <filename>

> The level data file describes the discrete levels and transitions of the outgoing (beam-like or target-like) reaction product. The file format is based on that of the PhotonEvaporationX.X data files (described in detail in the file $G4LEVELGAMMADATA/README-LevelGammaData). 
> 
> _Note: If the mass number of the outgoing nucleus lies outside of the A range specified in NuclearLevelData.cc by the AMIN[] and AMAX[] arrays, you will need to increase the range to include it so that your level data can be loaded._

Optional commands related to the outgoing reaction product:

    /BeamOut/XsectFile <filename>

> The differential cross section file is a text file containing a lab-frame scattering angle distribution. The format of the level scheme file:

>       <Minimum Theta [deg]> <Maximum Theta [deg]> <Theta bin width [deg]> 
>       <Channel 1 differential cross section [arbitrary units]>
>       <Channel 2 differential cross section [arbitrary units]>
>       <Channel 3 differential cross section [arbitrary units]>
>       ...

> If this file is present, the 2-body reaction will draw from this distribution to determine the scattering-angle for each event. The minimum and maximum scattering angles read from this file supersede values set with the `/BeamOut/ThetaMin` and `/BeamOut/ThetaMax` commands.

    /BeamOut/AngDistSigmaA <double> <unit>
    
    /BeamOut/AngDistSigmaB <double> <unit>

> Angular spreads of the lab-frame scattering angle distribution of the outgoing beam-like reaction products in the dispersive and nondispersive directions, respectively. The 2-body reaction kinematics draw from this scattering-angle distribution. This is the alternative to providing a lab-frame scattering-angle distribution with the `/BeamOut/XsectFile` command.

    /BeamOut/ThetaMin <double> <unit>
    
    /BeamOut/ThetaMax <double> <unit>

> Limits of the scattering angle distribution used in the 2-body reaction kinematics. If a scattering-angle distribution is supplied with the `/BeamOut/XsectFile` command, these values are superseded.

    /BeamOut/Source

> Simulate a stationary source using the in-beam simulation framework. The incoming beam particle decays into the outgoing beam particle. The decay of the incoming beam is handled by the  `G4RadioactiveDecay` process. A level scheme file must be used. The energy of the incoming beam must be set to zero, and the position of incoming beam must be set to the desired source position (see `./examples/sources/eu152/eu152.mac`,  `./examples/sources/co60`, and `./examples/sources/ho166`). _This is an independent approach to that described under **Source Simulations** below, in which gamma-rays are emitted as primary particles._ 

    /process/inactivate Reaction

> Turns off the Reaction process. This is required for simulations of stationary sources (with the /BeamOut/Source command) in which the `G4RadioactiveDecay` process manages gamma-ray emission. This can also be used for simulations of the beam passing through the target without reacting.

### Nuclear Alignment (Polarization) and Gamma-Ray Angular Distributions (see also ./examples/inbeam/angdist) ###

The alignment of the excited reaction product can be specified in in-beam simulations, leading to a net polarization of emitted gamma rays and a corresponding non-isotropic gamma-ray angular distribution. This functionality can be enabled by compiling with the `POL=1` flag (producing an executable named `UCGretina_Pol`).

The alignment of the reaction product must be specified (after the `/run/initialize` command) in terms of the population parameters P(m) of the magnetic substates, set using:

    /reaction/population <2m> <P(m)>

where m = -J, -J+1, ..., J-1, J and the P(m) sum to 1. 

It is also important to and include the commands (prior to the `/run/initialize` command):

    /PhysicsList/AngularCorrelations true
    /PhysicsList/SetGammaPolarization true

### Source Simulations (see also ./examples/eu152/eu152_gammas.mac) ###

Mandatory commands

    /Experiment/RunSource
    
    /Experiment/Source/Set <type>

> Currently implemented types: `eu152`, `cs137`, `co56`, `co60`, `ra226`, `am241`, `photopeaks`, `eu152_peaks`, `co56_peaks`, `ra226`, `au`, `white`, `simple`, `neutron`

> The simple source type emits gamma rays of a single energy (set with the `/Experiment/Source/setEnergy` command) 

> The white source type emits gamma rays in a uniform energy distribution set with the `/Experiment/Source/setWhiteLowE` and `/Experiment/Source/setWhiteHighE` commands. The multiplicity of the white source is set with the `/Experiment/Source/setMultiplicity` command.

> `eu152_peaks`, `co56_peaks`, and `photopeaks` sources produce selected gamma-rays from <SUP>152</SUP>Eu, <SUP>56</SUP>Co, and both, respectively. The gamma rays are emitted in equal quantities to facilitate determining photopeak efficiencies by fitting the photopeaks directly.

> _Note: The `eu152`, `co56`, `co60`, and `ra226` source types reproduce empirically-observed relative intensities. However, the total number of gamma rays simulated does not correspond to the number of decays of the source, because these simulated sources emit a single gamma-ray per event, while some gamma rays from these sources are emitted in cascades._

> The `neutron` source type emits neutrons instead of gamma rays, functioning in all other ways like the `simple` gamma-ray source.

Optional commands

    /Experiment/Source/setEnergy <double> <unit>

> Energy of the `simple` and `neutron` source types

    /Experiment/Source/setX <double> <unit>
    /Experiment/Source/setY <double> <unit>
    /Experiment/Source/setZ <double> <unit>

> Position of the source.

    /Experiment/Source/setR <double> <unit>

> Radius of the source disk. 

    /Experiment/Source/setDX <double> <unit>
    /Experiment/Source/setDY <double> <unit>

> Horizontal (nondispersive) and vertical (dispersive) widths of a rectangular source. These override the `/Experiment/Source/setR` command.

    /Experiment/Source/setSigmaX <double> <unit>
    /Experiment/Source/setSigmaY <double> <unit>

> Horizontal (nondispersive) and vertical (dispersive) sigma parameter of a Gaussian distribution of emission points. These override the `/Experiment/Source/setR` command. 

> Note: The setDX, setDY and setSigmaX, setSigmaY can be mixed (to give a flat dispersive and Gaussian nondispersive distribution of emission points, e.g.).

    /Experiment/Source/CollimationAngle <double> <unit>

> Angular spread of the collimated beam (about the collimation direction).

    /Experiment/Source/CollimationDirection <double> <double> <double>

> X, Y, and Z components of a vector specifying the direction of the collimated beam. 

   /Experiment/Source/ThetaFile <filename>

> File name for the theta distribution of the emitted particles. This is a text file with format: 

>       <Minimum theta [rad]> <Maximum theta [rad]> <theta bin width [rad]> 
>       <Bin 1 counts> 
>       <Bin 2 counts>
>       <Bin 3 counts>
>       ...

> If this command is present, the direction of each emitted particle is set by drawing randomly from the given theta distribution.

> (The above commands setting the position, radius, collimation, and theta distribution have no effect with the "white", "background", "bgwhite", and "muon" source types.)

    /Experiment/Source/setWhiteLowE  <double> <unit>
    
    /Experiment/Source/setWhiteHighE <double> <unit>

> Energy range of a flat distribution for the `white` and `bgwhite` source types (also work with `simple` and `neutron` source types, superseding setEnergy)

    /Experiment/Source/setMultiplicity <int>

> Multiplicity of the "white" and "bgwhite" source types

    /Target/sourceFrame <frame type>

> Optionally include the source frame (`/Target/Construct` command required).

> Currently implemented frames: `eu152_Z2707`, `cs137_E2879`, and `co56_2012`

### Background Simulations (see also ./examples/background) ###

Mandatory commands for running background simulations

    /Experiment/RunSource
    
    /Experiment/Source/Set < background || bgwhite || muon >

> The `background` source type emits several gamma rays from a solid spherical shell surrounding GRETINA. The `bgwhite` source type emits gamma rays from a uniform energy distribution from the solid shell. The `muon` source type emits 4.0 GeV muons vertically from above GRETINA.

Optional commands describing the spherical shell surrounding GRETINA from which background gamma-rays are emitted.

    /BackgroundSphere/Material <material name>

> Set material for the background sphere (default: `G4_Galactic`).

    /BackgroundSphere/R_min <double> <unit>
    
    /BackgroundSphere/R_max <double> <unit>

> Set the inner and outer radii of the background sphere (default: 3.0 m, 3.4 m). 

### LBL Scanning Table (see also ./examples/scan) ###

(`UCGretina_Scan` only)

 The geometries specified in `./GretinaGeometry/Scan0`, `./GretinaGeometry/Scan1`, `./GretinaGeometry/Scan2`, and `./GretinaGeometry/Scan3` orient the GRETINA module such that the corresponding crystal is centered on the slits.

Mandatory command for building the scanning table:

    /ScanningTable/CADModelPath <path>

> Path to the STL files (usually the cadModels directory of the UCGretina source tree).

Optional commands for the scanning table:

    /ScanningTable/Clover < left || right || both >

> Construct the clover detector(s).

    /ScanningTable/IncludeCartFrame

> Construct the scanning table frame and GRETINA mount.

    /ScanningTable/IncludeSlits

> Construct the slits.

    /ScanningTable/IncludeSlitMount

> Construct the slit mount.

    /ScanningTable/IncludeCollimator

> Construct the collimator.

    /ScanningTable/IncludeCollimatorInsert

> Construct the collimator insert.

    /ScanningTable/IncludeCollimatorMount

> Construct the collimator mount and x/y translation assemblies.

    /ScanningTable/IncludeCloverCart

> Construct the Clover cart.

    /ScanningTable/IncludeShields

> Construct the BGO anti-Compton shields.

    /ScanningTable/SetControllerX <double> <unit>
    
    /ScanningTable/SetControllerY <double> <unit>

> Set the horizontal positions of the source collimator relative to the central axis of the GRETINA module. (These positions correspond to those reported by the stepper motor controller. The controller x axis points opposite the geant4 x axis, and the controller y axis points along the geant4 z axis.) The horizontal position of the source should not be set using the usual source positioning commands (`/Experiment/Source/setX` and `/Experiment/Source/setZ`). 

    /ScanningTable/SetControllerZ <double> <unit>

> Set the vertical position of the slit assembly. (This position corresponds to that reported by the stepper motor controller. The controller z axis corresponds to the geant4 y axis.) 

    /ScanningTable/SetCloverZ <double> <unit>

> Set the vertical shift of the Clover detector assembly. (This position corresponds to the vertical position of the Clover assembly relative to its lowest position, determined by the brackets on the scanning table frame.)

    /ScanningTable/SetCollimatorRadius <double> <unit>

> Set the inner radius of the collimator insert.

### Liquid Hydrogen Target (see also ./examples/inbeam/fitLH) ###

(`UCGretina_LH` only)

Optional commands for setting LH target parameters (must precede /Target/Construct):

    /Target/Cell < thick || thin || empty || notarget >

> The "empty" type is the cell body with no window frames. (To construct an empty cell with frames, set "thick" or "thin" and set the target material to vacuum.) The "notarget" type constructs only the LH-target beam pipe without the target assembly.

    /Target/Bulge <double> <unit>

> Set the maximum target thickness (at the center) added to the target by the window bulge. (This increases the total target thickness along the beam axis by twice the bulge thickness.)

    /Target/Windows

> Include the Kapton cell windows (125 micron).

    /Target/Angle <double> <unit>

> Set the angle of tilt about the beam axis of the entire target assembly (30 degrees for GRETINA at the NSCL).

    /Target/SetDensity <double>

> Set target density in mg/cc.

    /Target/Material <material>

> Only `vacuum` or `G4_Galactic` are allowed with the LH target. (This is provided to enable source simulations with an empty cell with window frames installed.)

Mandatory command for building the LH target (`UCGretina_LH`):

    /Target/Construct

## Output ##

Output can be written in ASCII and/or "Mode 2" binary format. Raw tracking information can also be written to standard output.

### Position Resolution ###

The optional command

    /Output/PositionResolution <double> <unit>

sets the sigma parameter of a random Gaussian distribution folded into the simulated gamma-ray interaction positions (default = 0). 

### ASCII Output ###

The optional command

    /Output/Filename <filename>

activates ASCII output.

If energy was deposited in GRETINA in the event, S800 data and decomposed gamma-ray information are written: 

    S   <ATA>   <BTA>   <DTA>   <YTA>   <Event #>
    D   <# of decomposed gamma events>   <Event #>
    C   <Crystal ID>   <# of interaction points>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
    C   <Crystal ID>   <# of interaction points>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
    ...

The energy, emission position, emission direction, and the velocity of the projectile are written for each gamma ray emitted in every event:

    E   <# of emitted gamma rays>    <Full Energy>     <Event #>
        <Energy>   <X>   <Y>   <Z>   <theta>   <phi>   <beta>
        <Energy>   <X>   <Y>   <Z>   <theta>   <phi>   <beta>
        ...

`<Full Energy> = 1` if a single gamma ray is emitted and its full energy is deposited in a single crystal. `<Full Energy> = 0` if a single gamma ray is emitted and only part of its energy is deposited in any one crystal. `<Full Energy> = -1` otherwise.

The optional command

    /Output/DetectorsOnly

writes detected gamma-ray information only. Simulated S800 and emitted gamma-ray information is suppressed.

### Mode 2 Output ###

Mode 2 output from in-beam simulations contains S800 tracking events (GEB type 9), decomposed gamma-ray events (GEB type 1), and emitted gamma-ray (GEB type 11). (Source simulations do not produce S800 tracking events.) S800 tracking events are only written for events in which gamma rays are detected. Emitted gamma-ray events are written for every event. They include the energies, emission positions, and emission directions of all gamma rays emitted in an event.

Energies are expressed in keV, and positions are expressed in mm.

### Commands related to mode 2 output ###

    /Mode2/Filename <filename>

> This optional command activates Mode 2 output.

    /Mode2/crmatFile <filename>

> The file `crmat.LINUX` is provided. It specifies the rotation matrices and translation vectors from each crystal frame into the world coordinate system. UCGretina inverts this transformation to produce data in the crystal frames, as expected in Mode 2 data.

    /Mode2/crystalXforms

> Internal transformations from the world frame to the crystal frames are used. This is an alternative to specifying a crmat file. There are small deviations in crystal placement from the design positions. Therefore, using the standard crmat file introduces offsets in the crystal-frame hit patterns. Using the internal transformations does not.

    /Mode2/GretinaCoords

> When this command is present, the world coordinates are rotated Pi/2 about z to match the standard GRETINA coordinate system (x = down, z = beam). This is the coordinate system expected in Mode 2 data. This only affects the decomposed gamma-ray events. 

> Interaction points are expressed in the Geant4 coordinate system (y = up, z = beam) by default. 

    /Mode2/PackingRes <float> <unit>

> The PackingRes parameter determines the closest spacing of gamma-ray interaction points within each segment that GRETINA can resolve. Simulated interaction points are consolidated accordingly. Specify a zero position resolution to turn off consolidation.

    /Mode2/S800KE <float> <unit>

> Specifies the kinetic energy of the reaction product centered in the acceptance of the S800. This is needed for calculating DTA (dT/T) values of the S800 tracking events in the Mode2 output. 

    /Mode2/AllS800

> All S800 events are written, even in the absence of detected gamma rays.

    /Mode2/Print

> This command triggers printing of information on Mode 2 data to stdout. 

### Raw Tracking Information ###

    /GammaPrint/Track_Set

> Print gamma-ray tracking information to standard output.

    /IonPrint/Track_Set

> Print ion tracking information to standard output.

## Visualization ##

Run the macro file `vis/vis.mac` an interactive session:

    $ UCGretina
    
    Idle> /control/execute vis/vis.mac
    Idle> exit

This generates a VRML 2 file named `g4_XX.wrl` which can be viewed with a VRML viewer (like view3dscene, FreeWRL, or mayavi2). Macro files are included for visualizing GRETA, the LBL scanning table, and the liquid-hydrogen target setups. These scripts must be run in a directory with soft links (`aclust`, `euler`, `aslice`, `asolid`, `awalls`) to the appropriate GRETINA geometry files.

The macro file `./vis/trajectories.mac` illustrates how to add particle trajectories to visualizations.

Within mayavi2, the python scripts `./vis/mlab.animate.py` and `./vis/mlab.movie.py` can be run (File -> Run Python Script). The former animates the scene, and the latter saves the animation frames as a series of .png files which can be stitched together into an animated png or gif.
