# BrainGrowthModel
T.Tallinen_Codes_MA
Codes published by T. Tallinen: http://users.jyu.fi/~tutatall/codes.html (Brain Folder) Tallinen's articles:

On the growth and form of cortical convolutions, T. Tallinen et al, 2016 (https://www.nature.com/articles/nphys3632) (1)
Gyrification from constrained cortical expansion, T. Tallinen et al, 2014 (https://www.pnas.org/content/pnas/111/35/12667.full.pdf) (2)
Code has been modified with:

A normalisation of the mesh before starting the simulation, to set the dimensions of the brain longitudinal length of the mesh from -1 to +1. After running (folding & growing) but before generating output meshes, original dimensions are recovered by denormalising the mesh with the inverse of the normalisation factor. This is done to avoid conflicts in the Contact Processing section. This way the code works with any mesh of any size.

Timing in the model has been modified to also be in terms of gestational age (GA). In the original codes the parameter "t" ranges from 0 to 1 (meaning: 0 = 22 GA and 1 = 44 GA). At each step it is incremented by t =+ dt. As seen in the X. Wang's article (T. Tallinen also did that way): On early brain folding patterns using biomechanical growth modeling, 2019 (3): "t parametrizes time and has a non-linear relation to gestational age (GA) as t = 6.926 x 10^(-5)·GA^3 - 0.00665·GA^2 + 0.250·GA - 3:0189; where t [0,1] corresponds to GA [22 GA, adult]." *CHANGE EQUATION! Now model follows the Gompertz function proposed by X. Wang.

What has been done is to:

i) set manually the specific GA* of the specific input mesh. It should be different for different subjects;

ii) compute t* = t(GA*) with the equation proposed by X. Wang (and T. Tallinen);

iii) start running the simulation from this t* instead of t = 0. Run the simulation until t = 1. Simulations would run from GA* > 22 until GA == 44;

iv) (without further effect): in the simulation loop, at each step, GA is calculated again in terms of t. This is done with the inverse of the abovementioned equation. Output meshes and numerical results are generated every X GA, instead of every X steps.

Cortical thickness (H) Model: Defined in X. Wang's too (3). (*TO CHANGE!!! NOW IT IS DEFINED PATIENT-SPECIFIC+REGIONALLY). Instead of defining it constant throughout all the simulation, H is defined to increase steadily during the expanding process. H is defined as a linear function of time: H = Hi + B·t. where Hi is the initial cortical thickness and B (beta) determines the magnitude of thickening. Fetal brains of GA 22 have a typical cortical thickness 2.5 mm; the corresponding deformed cortical thickness in adulthood is approximately: 2.8 mm (t = 1). Then, B is defined to be B=0.3 in the simulations.

A function that generates simulated meshes in STL format.

Possibility of printing parameter values. Lines of code are commented by default, but if there is interest to see the values of different variables, lines can be uncommented and results are printed in a .dat file.

*TO COMPLETE