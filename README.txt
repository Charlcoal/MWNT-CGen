THIS CODE IS CURRENTLY NOT IMPLEMENTED IN AN ORGANIZED MANNER

Using MWNT-CGen requires editing the main function to call various parts of the code. MWNT-CGen currently only supports the Avogadro cjson file format, and will likely require exporting through Avogadro 2 to be used in other programs.

CITING:
If used in an academic paper, please cite github url.


IMPLEMENTATION:
(ALL UNITS IN nm UNLESS SPECIFIED)


||MWNTube(double innerDiam, double outerDiam, double length, double bondLength, double interWallDistance, bool fitOuterDiam (optional)):


|innerDiam: the target diameter of the inner-most layer of the tube. Will be exact if the tube has 1 layer and fitOuterDiam is set to false.

|outerDiam: the target diameter of the outer-most layer. Will be exact anytime fitOuterDiam is not explicitly set to false.

|bondLength: the distance between bonded carbon atoms. 0.1421 is the generally accepted value.

|interWallDistance: the difference bitween the radii of each layer, or the physical distance between the layers.

|fitOuterDiam: when true, the outerDiam will be exact. When false, the program will comprimise between inner and outer diameters, unless the interWallDistance is small enough for there to be 1 layer, in which case it will fit the inner diameter exactly.


|.generate(bool doubleBond): performes the computationally-heavy part of generating the nanotube. If doubleBond is true, then bonds running along the axis of the tube will be double-bonds. This is physically innacurate but may be necisary for some energy minimization programs.

|.addCatalyst(CrystalStructure structure, double radius, double density, double layerCatalystDist, double regularity): places hemispheres of an arbitrary crystal on the surface of the tube, WITHOUT BONDING. Hemispheres can overlap, but will not place atoms near already placed ones AS LONG AS .addCatalyst IS NOT CALLED MULTIPLE TIMES ON THE SAME TUBE.
    
    |structure: defines the crystal structure to be added. See the CrystalStructure implementation.
    
    |radius: the radius of hemispheres that contain the crystal.
    
    |density: the frequency of hemisphere generation, recomended value 0.1 -> 0.5
    
    |regularity: defined behavior for range [0, 1]. The higher the value, the less random hemisphere spacing will be. 0 is not entirely random, however, and still spreads out hemispheres somewhate uniformly.


|.output(ofstream outputFile): outputs the created nanotube (with any catalyst) to an ofstream. Requires ofstream declaration and opening before calling, as well as closing after calling. OFSTREAM OVERWRITES ANY EXISTING DATA IN THE FILE, AND .OUTPUT BREAKS IF CALLED MORE THAN ONCE ON THE SAME ofstream DURING THE SAME INSTANCE OF THE PROGRAM (running the program as a whole twice will overwrite the original).



||CrystalStructure(double alpha, double beta, double gamma, double a, double b, double c, vector<Vec3d> fractionalAtomPos, vector<int> atomElementNumbers):


|alpha, beta, and gamma: triclinic crystal cell angles (in degrees)

|a, b, and c: triclinic side lengths

|fractionalAtomPos: the position of crystal elements in the triclinic coordinate system. Vec3d is constructed as {x, y, z}. For best results, ensure that the face/edge/corner atoms are never repeated (eg. don't specify both {0, 0, 0} and {1, 0, 0}), and are placed as close to other atoms in the cell as possible (not required).

|atomElementNumbers: the atomic numbers of the corisponding atom in the fractionAtomPos (eg. 1 for hydrogen and 8 for oxygen).
