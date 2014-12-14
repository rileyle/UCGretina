# UCGretina v2.9.1 #

## Compile and install ##

Install version 4.9.6.p03 of the Geant4 libraries from
http://geant4.web.cern.ch/geant4/support/download.shtml. You will need
the data files for low energy electromagnetic processes. 

_Important note: UCGretina is not yet compatible with the latest
version (4.10.x) of the Geant4 libraries._

Set up your environment (consider adding this to your .bashrc):

    $ source <G4INSTALL>/share/Geant4-9.6.3/geant4make/geant4make.sh

Compile:

    $ make

To use the liquid hydrogen target:

    $ make LHTARGET=1

(produces the binary UCGretina_LH)

The executables UCGretina and UCGretina_LH are automatically installed
in

    $G4WORKDIR/bin/$G4SYSTEM

(which is added to your path when you source geant4make.sh)

## Examples ##

Several examples are in the examples subdirectory, including
illustrations of fitting simulations to measured source and in-beam
spectra. In the examples, the GrROOT package is used for sorting the
mode 2 output from UCGretina.

## Selected Macro File Commands ##

### Geometry ###

Optional commands for setting target parameters:

    /Target/Material <material name>

    /Target/X_length <double> <unit>

    /Target/Y_length <double> <unit>

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

Mandatory command for building the GRETA chamber (spherical chamber
with 36 cm OD and 1.5 mm wall thickness, and a beam pipe with 4 cm OD
and 1.5 mm wall thickness):

    /GretaChamber/Construct

Optional commands for setting GRETA chamber geometry:

    /GretaChamber/R_min <double> <unit>

    /GretaChamber/R_max <double> <unit>

Optional commands for including Gretina-related dead material:

    /Gretina/detector/enableCapsules

    /Gretina/detector/enableCryostats

    /Gretina/Shell < full || north || south || Greta || GretaLH || Greta_North || Greta_South || GretaLH_North || GretaLH_South >

Mandatory command after setting any GRETINA parameters:

    /Gretina/update

### LH Target Geometry ###

(UCGretina_LH only)

Optional commands for setting LH target parameters (must precede
/Target/Construct):

    /Target/Cell < thick || thin || empty || notarget >

> The "empty" type is the cell body with no window frames. (To
  construct an empty cell with frames, set "thick" or "thin" and set
  the target material to vacuum.) The "notarget" type constructs only
  the LH-target beam pipe without the target assembly.

    /Target/Bulge <double> <unit>

> Set the maximum target thickness (at the center) added to the target
  by the window bulge. (This increases the total target thickness
  along the beam axis by twice the bulge thickness.)

    /Target/Windows

> Include the Kapton cell windows (125 micron).

    /Target/Angle <double> <unit>

> Set the angle of tilt about the beam axis of the entire target
  assembly (30 degrees for Gretina at the NSCL).

    /Target/SetDensity <double>

> Set target density in mg/cc.

    /Target/Material <material>

> Only "vacuum" or "G4_Galactic" are allowed with the LH target. (This
  is provided to enable source simulations with an empty cell with
  window frames installed.)

Mandatory command for building the LH target (UCGretina_LH):

    /Target/Construct

### In-beam Simulations ###

Commands related to the incoming beam:

    /BeamIn/A <A> /BeamIn/Z <Z>

> Mass number and atomic number

    /BeamIn/KEu <double> <unit>

> Kinetic energy per nucleon of the incoming beam

    /BeamIn/Dpp <double>

> Momentum acceptance (dp/p) for the incoming beam. (Set this to 0 if
> you are providing the dta spectrum of the incoming beam using the
> /BeamIn/dtaFile command.)

    /BeamIn/dtaFile

> file name for the dta spectrum of the incoming beam. This is a text
> file with format: 

>       <Minimum DTA [%]> <Maximum DTA [%]> <DTA bin width [%]> 
>       <Channel 1 counts> 
>       <Channel 2 counts>
>       <Channel 3 counts>
>       ...

> If this command is present, the incoming beam energy is set by
> drawing randomly from the dta distribution centered on the beam
> energy specified by the /BeamIn/KEu command. (The dta spectrum
> specifies the momentum acceptance of the incoming beam, so the
> momentum acceptance parameter should be set to zero: /BeamIn/Dpp 0.)

    /BeamIn/Focus/X <double> <unit> /BeamIn/Focus/Y <double> <unit>

> Offsets of the emission point of the incoming beam. (Z defaults to
> -50 cm. If you change this, make sure it is upstream of the target!)

    /BeamIn/Focus/DX <double> <unit> /BeamIn/Focus/DY <double> <unit>

> Horizontal and vertical widths of the beam spot at the emission
> point (not on target)

    /BeamIn/Focus/Ata0 <double> <unit> /BeamIn/Focus/Bta0 <double> <unit>

> Direction of the incoming beam (dispersive and nondispersive angles,
> respectively)

    /BeamIn/Focus/maxAta <double> <unit> 

    /BeamIn/Focus/maxBta <double> <unit>

> Angular divergences of the incoming beam in the dispersive and
> nondispersive directions, respectively.

Commands related to the outgoing reaction product:

    /BeamOut/DA <int> 

    /BeamOut/DZ <int>

> Changes in mass number and atomic number of the reaction. The
> incoming beam has (A,Z) and the outgoing reaction product has (A+DA,
> Z+DZ)

    /BeamOut/TargetA <int>

    /BeamOut/TargetZ <int>

> Mass number and charge number of the target nucleus.

    /BeamOut/ProjectileExcitation <double> <unit>

> Excitation energy of the beam-like reaction product. (The
> target-like reaction product is not excited.) If a level scheme file
> is used, the energy parameter is ignored. 

    /BeamOut/TargetExcitation <double> <unit>

> Excitation energy of the target-like reaction product. (The
> beam-like reaction product is not excited.) The energy parameter is
> superseded by a level scheme file if a /BeamOut/LevelSchemeFile
> command is present.

    /BeamOut/XsectFile <filename>

> The differential cross section file is a text file containing a
> lab-frame scattering angle distribution. 
> The format of the level scheme file:

>       <Minimum Theta [deg]> <Maximum Theta [deg]> <Theta bin width [deg]> 
>       <Channel 1 differential cross section [arbitrary units]>
>       <Channel 2 differential cross section [arbitrary units]>
>       <Channel 3 differential cross section [arbitrary units]>
>       ...

> If this file is presnet, the 2-body reaction will draw from this
> distribution to determine the scattering-angle for each event. The
> minimum and maximum scattering angles read from this file supersede
> values set with the /BeamOut/ThetaMin and /BeamOut/ThetaMax
> commands.

    /BeamOut/AngDistSigmaA <double> <unit>

    /BeamOut/AngDistSigmaB <double> <unit>

> Angular spreads of the lab-frame scattering angle distribution of
> the outgoing beam-like reaction products in the dispersive 
> and nondispersive directions, respectively. The 2-body reaction
> kinematics draw from this scattering-angle distribution. This is
> the alternative to providing a lab-frame scattering-angle distribution 
> with the /BeamOut/XsectFile command.

    /BeamOut/ThetaMin <double> <unit>

    /BeamOut/ThetaMax <double> <unit>

> Limits of the scattering angle distribution used in the 2-body
> reaction kinematics. If a scattering-angle distribution is supplied
> with the /BeamOut/XsectFile command, these values are superseded.

    /BeamOut/seta0 <double> 

    /BeamOut/seta2 <double> 

    /BeamOut/seta4 <double>

> Angular distribution coefficients for the emitted gamma rays.
> These commands are superseded by a level scheme file if a
> /BeamOut/LevelSchemeFile command is present.

    /BeamOut/tau <double> <unit>

> Mean lifetime of the excitation. This command is superseded by a
> level scheme file if a /BeamOut/LevelSchemeFile command is present.

    /BeamOut/LevelSchemeFile <filename>

> The level scheme file describes the portion of the level scheme to
>  be simulated. The format of the level scheme file:

>       <Level energy [keV]> <Number of gamma-decay Branches> <Lifetime [ps]> <Relative direct population> 
>       <BR 1> <Final-state energy [keV]> <A0> <A2> <A4> 
>       <BR 2> <Final-state energy [keV]> <A0> <A2> <A4> 
>       ...  
>       <Level energy [keV]> <Number of gamma-decay Branches> <Lifetime [ps]> <Relative direct population> 
>       <BR 1> <Final-state energy [keV]> <A0> <A2> <A4>
>       <BR 2> <Final-state energy [keV]> <A0> <A2> <A4>
>       ...

> where <BR N> are branching ratios, <A0>, <A2>, <A4> are gamma-ray
> angular distribution coefficients. The <Relative direct population>
> parameters determine how often the level is populated directly by
> the reaction.

> NOTE: This command supersedes values set with the
> /BeamOut/ProjectileExcitation, /BeamOut/TargetExcitation,
> /BeamOut/seta0, /BeamOut/seta2, /BeamOut/seta4, and /BeamOut/tau
> commands described above. 

    /BeamOut/Source

> Simulate a stationary source using the in-beam simulation framework.
> The incoming beam particle decays into the outgoing beam particle. A
> level scheme file must be used. The energy of the incoming beam must
> be set to zero, and the position of incoming beam must be set to the
> desired source position. _This is an independent approach to that
> described under **Source Simulations** below, in which gamma-rays are 
> emitted as primary particles._ 

### Source Simulations (see also ./examples/eu152) ###

Mandatory commands

    /Experiment/RunSource

    /Experiment/Source/Set <type>

> Currently implemented types: eu152, cs137, co56, co60,
> photopeaks, eu152_peaks, co56_peaks, au, white, simple

> The simple source type emits gamma rays of a single energy (set with
> the /Experiment/Source/setEnergy command) 

> The white source type emits gamma rays in a uniform energy
> distribution between 100 keV and 10 MeV. These limits can be set
> with the /Experiment/Source/setWhiteLowE and
> /Experiment/Source/setWhiteHighE commands.

> eu152_peaks, co56_peaks, and photopeaks sources produce selected
> gamma-rays from 152Eu, 56Co, and both, respectively. The gamma rays
> are emitted in equal quantities to facilitiate determining photopeak
> efficiencies by fitting the photopeaks directly.

Optional commands

    /Experiment/Source/setEnergy <double> <unit>

> Set the energy of the "simple" source type

    /Experiment/Source/setWhiteLowE  <double> <unit>

    /Experiment/Source/setWhiteHighE <double> <unit>

> Set the limits of the energy distribution of the "white" source type

    /Target/sourceFrame <frame type>

> Optionally include the source frame (/Target/Construct command required).

> Currently implemented frames: eu152_Z2707, cs137_E2879, and co56_2012

### Background Simulations (see also ./examples/background) ###

Mandatory commands for running background simulations

    /Experiment/RunSource

    /Experiment/Source/Set <background || bgwhite || muon>

> The background source type emits several gamma rays from a solid
> spherical shell surrounding GRETINA. The bgwhite source type emits
> gamma rays from a uniform energy distribution from the solid
> shell. The muon source type emits 4.0 GeV muons vertically from
> above GRETINA. 

Optional commands describing the spherical shell surrounding GRETINA from which background gamma-rays are emitted.

    /BackgroundSphere/Material <material name>

> Set material for the background sphere (default: G4_Galactic).

    /BackgroundSphere/R_min <double> <unit>

    /BackgroundSphere/R_max <double> <unit>

> Set the inner and outer radii of the background sphere 
> (default: 3.0 m, 3.4 m). 

## Tracking ##

    /GammaPrint/Track_Set

> Print gamma-ray tracking information to standard output.

    /IonPrint/Track_Set

> Print ion tracking information to standard output.

## Mode 2 Output ##

Mode 2 output from in-beam simulations contains S800 tracking events
(GEB type 9), decomposed gamma-ray events (GEB type 1), and emitted
gamma-ray (GEB type 11). (Source simulations do not produce S800
tracking events.) S800 tracking events are only written for events in
which gamma rays are detected. Emitted gamma-ray events are written
for every event. They include the energies, emission positions, and
emission directions of all gamma rays emitted in an event.

Energies are expressed in keV, and positions are expressed in mm.

### Commands related to mode 2 output ###

    /Mode2/Filename <filename>

> This optional command activates Mode 2 output.

    /Mode2/crmatFile <filename>

> The file "crmat.LINUX" is provided. It specifies the rotation
> matrices and translation vectors from each crystal frame into the
> world coordinate system. UCGretina inverts this transformation to
> produce data in the crystal frames, as expected in Mode 2 data.

    /Mode2/crystalXforms

> Internal transformations from the world frame to the crystal 
> frames are used. This is an alternative to specifying a crmat
> file. There are small deviations in crystal placement from the
> design positions. Therefore, using the standard crmat file
> introduces offsets in the crystal-frame hit patterns. Using the
> internal transormations does not.

    /Mode2/GretinaCoords

> When this command is present, the world coordinates are rotated Pi/2
> about z to match the standard GRETINA coordinate system (x = down, z
> = beam). This is the coordinate system expected in Mode 2 data. This
> only affects the decomposed gamma-ray events. 

> Interaction points are expressed in the Geant4 coordinate system 
> (y = up, z = beam) by default. 

    /Mode2/PackingRes <float> <unit>

> The PackingRes parameter determines the closest spacing of
> gamma-ray interaction points within each segment that GRETINA can
> resolve. Simulated interaction points are consolidated
> accordingly. Specify a zero position resolution to turn off
> consolidation.

    /Mode2/S800KE <float> <unit>

> Specifies the kinetic energy of the reaction product centered in the
> acceptance of the S800. This is needed for calculating DTA
> (dT/T) values of the S800 tracking events in the Mode2 output. 

    /Mode2/Print

> This command triggers printing of information on Mode 2 data to
> stdout. 

### ASCII Output ###

The optional command

    /Output/Filename <filename>

activates ASCII output.

If energy was deposited in GRETINA in the event, S800 data and
decomposed gamma-ray information are written: 

    S   <ATA>   <BTA>   <DTA>   <YTA>   <Event #>
    D   <# of decomposed gamma events>   <Event #>
    C   <Crystal ID>   <# of interaction points>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
    C   <Crystal ID>   <# of interaction points>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
        <Segment ID>   <Energy>   <X>   <Y>   <Z>
    ...

The energy, emission position, emission direction, and the
velocity of the projectile are written for each gamma ray emitted in
every event:

    E   <# of emitted gamma rays>    <Full Energy>     <Event #>
        <Energy>   <X>   <Y>   <Z>   <theta>   <phi>   <beta>
        <Energy>   <X>   <Y>   <Z>   <theta>   <phi>   <beta>
        ...

<Full Energy> = 1 if a single gamma ray is emitted and its
full energy is deposited in a single crystal. <Full Energy> = 0 if a
single gamma ray is emitted and only part of its energy is deposited
in any one crystal. <Full Energy> = -1 otherwise.

## Visualize the array ##

Run the macro file "vis/vis.mac" an interactive session:

    $ UCGretina

    Idle> /control/execute vis/vis.mac
    Idle> exit

This generates a VRML 2 file named g4_XX.wrl which can be viewed
with a VRML viewer (like mayavi2).

Within mayavi2, the python scripts ./vis/mlab.animate.py and
./vis/mlab.movie.py can be run (File -> Run Python Script). The
former animates the scene, and the latter saves the animation frames
as a series of .png files which can be stitched together into an
animated png or gif.

