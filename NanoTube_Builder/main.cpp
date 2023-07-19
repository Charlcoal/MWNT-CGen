//MORE INFO & LATEST VERSION AT:
//https://github.com/Charlcoal/MWNT-CGen

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

class Vec3d {
public:
    double x = 0;
    double y = 0;
    double z = 0;
public:
    double magnitude();
    double magSquared();
};

double Vec3d::magSquared() {
    return (x * x) + (y * y) + (z * z);
}

double Vec3d::magnitude() {
    return sqrt(magSquared());
}

class Vec2d {
public:
    double x = 0;
    double y = 0;
public:
    double magnitude();
    double magSquared();
};

double Vec2d::magSquared() {
    return (x * x) + (y * y);
}

double Vec2d::magnitude() {
    return sqrt(magSquared());
}

double floatMod(double num, double modifier) { //WARNING: SLOW (try not to input floats with a high quotient)
    while(num < 0) {
        num += modifier;
    }

    while(num >= modifier) {
        num -= modifier;
    }
    return num;
}

Vec2d operator+(Vec2d const& a, Vec2d const& b) {
    return Vec2d{a.x + b.x, a.y + b.y};

}

Vec2d operator-(Vec2d const& a, Vec2d const& b) {
    return Vec2d{a.x - b.x, a.y - b.y};

}

Vec3d operator-(Vec3d const& a, Vec3d const& b) {
    return Vec3d{a.x - b.x, a.y - b.y, a.z - b.z};

}

Vec3d operator+(Vec3d const& a, Vec3d const& b) {
    return Vec3d{a.x + b.x, a.y + b.y, a.z + b.z};

}

double layerDiamCalculate(int n, int m, double bondLength) {
    return bondLength * sqrt(3 * ((n*n) + (m*m) + (n*m))) / M_PI;
}

class MWNTube;
//UNITS IN nm UNLESS SPECIFIED
class CNTLayer {
private:
    double bondLength;
    vector<Vec3d> atomPos;
    vector<pair<int, int>> bonds;
    vector<int> bondOrder;

    friend MWNTube;
public:
    double length;
    int atomLength;
    int n;
    int m;
    double chiralitySlope;
    double diameter;

private:
    int deriveAtomLength();
    void generateHex(vector<vector<pair<int, int>>>& dualGraphenePos);
    void addHexToAtoms(vector<vector<bool>>& yAtomPos, vector<vector<bool>>& deltaAtomPos, const vector<vector<pair<int, int>>>& dualGraphenePos, int atomArrayWidth);
    void determineBonds(const vector<vector<bool>>& yAtomPos, const vector<vector<bool>>& deltaAtomPos, int atomArrayWidth, int atomArrayHeight, vector<vector<vector<bool>>>& yBond);
    void placeBonds(bool doubleBond, int atomArrayHeight, int atomArrayWidth, const vector<vector<vector<bool>>>& yBond, vector<vector<bool>>& deltaAtomPos, vector<Vec2d>& atomFlatPos);
    void wrap(vector<Vec2d>& atomFlatPos);
public:
    void generate(bool doubleBond);
    CNTLayer(int n, int m, double length, double bondLength);
    void outPut();
};


int CNTLayer::deriveAtomLength() {
    return length / (1.5 * sqrt(1 + (chiralitySlope * chiralitySlope)) * bondLength);
}

CNTLayer::CNTLayer(int N, int M, double Length, double BondLength) {
    bondLength = BondLength;
    if(N < M) {
        n = M;
        m = N;
    } else {
        n = N;
        m = M;
    }
    length = Length;
    chiralitySlope = sqrt(3)*m/(2*n + m);
    diameter = layerDiamCalculate(N, M, bondLength);
    atomLength = deriveAtomLength();
}


vector<pair<int, int>> generateLoop(int n, int m) {
    vector<pair<int, int>> loop;
    int s = n + m;
    int i = 0;
    int j = 0;
    while(i + j < s) {
        loop.push_back({i, j});
        if(i == 0 || j/double(i) >= m/double(n)) {
            i++;
        } else {
            j++;
        }
    }

    return loop;
}

void CNTLayer::generateHex(vector<vector<pair<int, int> > > &dualGraphenePos) {
    dualGraphenePos.push_back(generateLoop(n, m));

    for(int i = 1; i < atomLength + 1; i++) {
        dualGraphenePos.push_back({});
        for(int j = 0; j < n + m; j++) {
            dualGraphenePos[i].push_back({dualGraphenePos[i - 1][j].first - 1, dualGraphenePos[i - 1][j].second + 1});
        }
    }
}

void CNTLayer::addHexToAtoms(vector<vector<bool> > &yAtomPos, vector<vector<bool> > &deltaAtomPos, const vector<vector<pair<int, int> > > &dualGraphenePos, int atomArrayWidth) {
    for(int i = 0; i < dualGraphenePos.size(); i++) {
        for(int j = 0; j < dualGraphenePos[i].size(); j++) {
            pair<int, int> center = {dualGraphenePos[i][j].second + 1, dualGraphenePos[i][j].first + dualGraphenePos[i][j].second};

            vector<pair<int, int>> offsets = {
                {0, 0},
                {1, 0},
                {1, 1}
            };

            for(int i = 0; i < offsets.size(); i++) {
                pair<int, int> point = {(center.first + offsets[i].first - (center.second + offsets[i].second >= atomArrayWidth ? m : 0)), (center.second + offsets[i].second) % atomArrayWidth};

                deltaAtomPos[point.first][point.second] = true;
            }

            for(int i = 0; i < offsets.size(); i++) {
                pair<int, int> point = {(center.first + offsets[i].second - (center.second + offsets[i].first >= atomArrayWidth ? m : 0)), (center.second + offsets[i].first) % atomArrayWidth};

                yAtomPos[point.first][point.second] = true;
            }
        }
    }
}

void CNTLayer::determineBonds(const vector<vector<bool> > &yAtomPos, const vector<vector<bool> > &deltaAtomPos, int atomArrayWidth, int atomArrayHeight, vector<vector<vector<bool> > > &yBond) {
    for(int i = 0; i < atomArrayHeight; i++) {
        for(int j = 0; j < atomArrayWidth; j++) {
            if(yAtomPos[i][j]) {
                yBond[i][j][0] = deltaAtomPos[i][j]; //top right segment
                yBond[i][j][1] = ((i + 1) < atomArrayHeight) && deltaAtomPos[i + 1][j]; //bottom segment
                if(j >= 1 || i + m < atomArrayHeight) {
                    yBond[i][j][2] = deltaAtomPos[i + (j < 1 ? m : 0)][(j - 1 + atomArrayWidth)%atomArrayWidth]; //top left segment
                }
            }
        }
    }
}

void CNTLayer::placeBonds(bool doubleBond, int atomArrayHeight, int atomArrayWidth, const vector<vector<vector<bool> > > &yBond, vector<vector<bool> > &deltaAtomPos, vector<Vec2d> &atomFlatPos) {
    vector<vector<int>> deltaAtomIndex(atomArrayHeight, vector<int>(atomArrayWidth));

    for(int i = 0; i < atomArrayHeight; i++) {
        for(int j = 0; j < atomArrayWidth; j++) {
            int count = 0;
            if(yBond[i][j][0]) {
                atomFlatPos.push_back(Vec2d{(j - (i * 0.5)) * sqrt(3) * bondLength, 1.5 * i * bondLength});
                count = 2;
                bondOrder.push_back(1);

                if(deltaAtomPos[i][j]) {
                    atomFlatPos.push_back(Vec2d{(j - (i * 0.5) + 0.5) * sqrt(3) * bondLength, (1.5 * i - 0.5) * bondLength});
                    deltaAtomPos[i][j] = false;
                    deltaAtomIndex[i][j] = atomFlatPos.size() - 1;

                    bonds.push_back({atomFlatPos.size() - count, atomFlatPos.size() - 1});
                } else {
                    count--;
                    bonds.push_back({atomFlatPos.size() - count, deltaAtomIndex[i][j]});
                }
            }

            if(yBond[i][j][1] && i + 1 < atomArrayHeight) {
                if(count == 0) {
                    atomFlatPos.push_back(Vec2d{(j - (i * 0.5)) * sqrt(3) * bondLength, 1.5 * i * bondLength});
                    count = 2;
                } else {
                    count ++;
                }

                bondOrder.push_back(1 + (doubleBond ? 1: 0)); //2 if doubleBonding is happening

                if(deltaAtomPos[i + 1][j]) {
                    atomFlatPos.push_back(Vec2d{(j - (i * 0.5)) * sqrt(3) * bondLength, (1.5 * i + 1) * bondLength});
                    deltaAtomPos[i + 1][j] = false;
                    deltaAtomIndex[i + 1][j] = atomFlatPos.size() - 1;

                    bonds.push_back({atomFlatPos.size() - count, atomFlatPos.size() - 1});
                } else {
                    count--;
                    bonds.push_back({atomFlatPos.size() - count, deltaAtomIndex[i + 1][j]});
                }
            }

            if(yBond[i][j][2]) {
                if(count == 0) {
                    atomFlatPos.push_back(Vec2d{(j - (i * 0.5)) * sqrt(3) * bondLength, 1.5 * i * bondLength});
                    count = 2;
                } else {
                    count += 1;
                }
                bondOrder.push_back(1);

                int jNeo = (j - 1 + atomArrayWidth) % atomArrayWidth;
                int iNeo = i + (j < 1 ? m : 0);
                if(deltaAtomPos[iNeo][jNeo]) {
                    atomFlatPos.push_back(Vec2d{(jNeo - (iNeo * 0.5) + 0.5) * sqrt(3) * bondLength, (1.5 * iNeo - 0.5) * bondLength});
                    deltaAtomPos[iNeo][jNeo] = false;
                    deltaAtomIndex[iNeo][jNeo] = atomFlatPos.size() - 1;

                    bonds.push_back({atomFlatPos.size() - count, atomFlatPos.size() - 1});
                } else {
                    count--;
                    bonds.push_back({atomFlatPos.size() - count, deltaAtomIndex[iNeo][jNeo]});
                }
            }
        }
    }
}

void CNTLayer::wrap(vector<Vec2d> &atomFlatPos) {
    Vec2d c{(n + m/2.0) * sqrt(3), m * 1.5};
    double mag = sqrt(c.x * c.x + (c.y * c.y));
    c = {c.x/mag, c.y/mag};
    Vec2d cPerp{-c.y, c.x};

    //Vec2ds are the columns
    double det = (c.x * cPerp.y) - (c.y * cPerp.x);
    pair<Vec2d, Vec2d> transMatrix = {{cPerp.y/det, -1 * c.y/det}, {-1 * cPerp.x/det, c.x/det}};

    for(int i = 0; i < atomFlatPos.size(); i++) { //transforms so that c is the unitvector for x
        atomFlatPos[i] = {transMatrix.first.x * atomFlatPos[i].x + transMatrix.second.x * atomFlatPos[i].y, transMatrix.first.y * atomFlatPos[i].x + transMatrix.second.y * atomFlatPos[i].y};

        atomFlatPos[i].x = 2 * atomFlatPos[i].x/(diameter); //turns c value into rads angle in cylinder

        atomPos.push_back({diameter * cos(atomFlatPos[i].x)/2, diameter * sin(atomFlatPos[i].x)/2, atomFlatPos[i].y});
    }
}

void CNTLayer::generate(bool doubleBond = false) {
    if(n < 1 || m < 0 || atomLength < 0) {
        return;
    }

    vector<vector<pair<int, int>>> dualGraphenePos; //dual center (triangle) positions in graphene (n, m)

    generateHex(dualGraphenePos);

    int atomArrayHeight = atomLength + m + 3;
    int atomArrayWidth = n + m;

    vector<vector<bool>> yAtomPos(atomArrayHeight, vector<bool>(atomArrayWidth, false));
    vector<vector<bool>> deltaAtomPos(atomArrayHeight, vector<bool>(atomArrayWidth, false));

    addHexToAtoms(yAtomPos, deltaAtomPos, dualGraphenePos, atomArrayWidth);

    vector<Vec2d> atomFlatPos;
    vector<vector<vector<bool>>> yBond(atomArrayHeight, vector<vector<bool>>(atomArrayWidth, vector<bool>(3, false)));

    determineBonds(yAtomPos, deltaAtomPos, atomArrayWidth, atomArrayHeight, yBond);
    placeBonds(doubleBond, atomArrayHeight, atomArrayWidth, yBond, deltaAtomPos, atomFlatPos);
    wrap(atomFlatPos);
}




Vec3d triclinicToEuclidean(Vec3d vec, double alpha, double beta, double gamma, double a, double b, double c) {
    double cx = cos(beta);
    double cy = (cos(alpha) - (cos(beta) * cos(gamma))) / sin(gamma);
    double z = sqrt(1 - (cx * cx) - (cy * cy));

    cx *= vec.z * c;
    cy *= vec.z * c;
    z *= vec.z * c;

    double x = vec.x * a + (vec.y * b * cos(gamma)) + cx;
    double y = vec.y * b * sin(gamma) + cy;

    return {x, y, z};
}

class CrystalStructure {
public:
    double alpha;
    double beta;
    double gamma;
    double a;
    double b;
    double c;
    vector<Vec3d> fractionalAtomPos;
    vector<int> atomElementNumbers;
public:
    Vec3d atomPos(int atomIndex);
};

Vec3d CrystalStructure::atomPos(int atomIndex) {
    return triclinicToEuclidean(fractionalAtomPos[atomIndex], alpha, beta, gamma, a, b, c);
}


class MWNTube {
private:
    vector<CNTLayer> layers;
    vector<pair<int, int>> layerParams;
    vector<Vec3d> catalystAtomPos;
    vector<int> catalystElementNumbers;
private:
    pair<int, int> generateLayerParam(pair<int, int> guess, double targetDiam, double bondLength);
public:
    MWNTube(double innerDiam, double outerDiam, double Length, double bondLength, double interLayerDistance, bool fitOuterDiam);
    void generate(bool doubleBond);
    int sp2Atoms();
    int numAtoms();
    bool hasBonds();
    void output(ofstream& outputFile);
    void addCatalyst(CrystalStructure structure, double radius, double density, double space, double regularity);
};


pair<int, int> MWNTube::generateLayerParam(pair<int, int> guess, double targetDiam, double bondLength) {
    vector<pair<int, int>> offsets{{-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, 1}, {1, -1}, {-1, 2}, {1, -2}, {-2, 1}, {2, -1}, {2, -3}, {-2, 3}, {3, -2}, {-3, 2}};
    bool repeat = true;
    pair<int, int> runningGuess = guess;
    while(repeat) {
        repeat = false;
        double stdError = abs(layerDiamCalculate(guess.first, guess.second, bondLength) - targetDiam);
        for(int i = 0; i < offsets.size(); i++) {
            if(guess.second + offsets[i].second < 0 ||
               guess.second + offsets[i].second > guess.first + offsets[i].first ||
               guess.first + guess.second + offsets[i].first + offsets[i].second < 3) {
                continue;
            }

            double error = abs(layerDiamCalculate(guess.first + offsets[i].first, guess.second + offsets[i].second, bondLength) - targetDiam);
            if(error < stdError) {
                stdError = error;
                repeat = true;
                runningGuess = {guess.first + offsets[i].first, guess.second + offsets[i].second};
            }
        }
        guess = runningGuess;
    }
    return guess;
}

MWNTube::MWNTube(double innerDiam, double outerDiam, double length, double bondLength, double interWallDistance, bool fitOuterDiam = true) {
    double diameterDifference = 2 * interWallDistance;

    if(fitOuterDiam) {
        vector<pair<int, int>> reversedLayerParams;
        reversedLayerParams.push_back(generateLayerParam({outerDiam * 2.4/bondLength, 1.2 * outerDiam/bondLength}, outerDiam, bondLength));
        double runningDiam = outerDiam - diameterDifference;

        while(runningDiam > innerDiam && (1.45868 * bondLength) - (runningDiam - diameterDifference) < (0.2 * diameterDifference)) {
            reversedLayerParams.push_back(generateLayerParam({runningDiam * 2.4/bondLength, 1.2 * runningDiam/bondLength}, runningDiam, bondLength));
            runningDiam-= diameterDifference;
        }

        for(int i = reversedLayerParams.size() - 1; i >= 0; i--) {
            layerParams.push_back(reversedLayerParams[i]);
        }
    } else {
        if((outerDiam - innerDiam / diameterDifference) > 1) {
            if(floatMod(outerDiam - innerDiam, diameterDifference) > diameterDifference/2) {
                innerDiam += (diameterDifference - floatMod(outerDiam - innerDiam, diameterDifference))/2;
            } else {
                innerDiam -= floatMod(outerDiam - innerDiam, diameterDifference)/2;
            }
        }

        layerParams.push_back(generateLayerParam({innerDiam * 2.4/bondLength, 1.2 * innerDiam/bondLength}, innerDiam, bondLength));
        double runningDiam = innerDiam + diameterDifference;

        while(runningDiam < outerDiam) {
            layerParams.push_back(generateLayerParam({runningDiam * 2.4/bondLength, 1.2 * runningDiam/bondLength}, runningDiam, bondLength));
            runningDiam += diameterDifference;
        }
    }


    for(int i = 0; i < layerParams.size(); i++) {
        layers.push_back(CNTLayer{layerParams[i].first, layerParams[i].second, length, bondLength});
    }
}

void MWNTube::generate(bool doubleBond = false) {
    for(int i = 0; i < layers.size(); i++) {
        layers[i].generate(doubleBond);
    }
}

int MWNTube::sp2Atoms() {
    int total = 0;
    for(int i = 0; i < layerParams.size(); i++) {
        total += layerParams[i].first;
        total += layerParams[i].second;
    }
    return total;
}

bool MWNTube::hasBonds() {
    for(int i = 0; i < layers.size(); i++) {
        if(layers[i].bonds.size() > 0) {
            return true;
        }
    }
    return false;
}

int MWNTube::numAtoms() {
    int count = 0;

    for(int i = 0; i < layers.size(); i++) {
        count += layers[i].atomPos.size();
    }

    return count;
}

void MWNTube::output(ofstream& outputFile) {
    outputFile << "{\n  \"chemicalJson\": 1, \n  \"atoms\": {\n    \"elements\": {\n      \"number\": [";

    if(numAtoms() > 0) {
        for(int i = 0; i < numAtoms(); i++) {
            outputFile << "  6" << ((i < numAtoms() - 1) ? ",": "");
        }

        for(int i = 0; i < catalystElementNumbers.size(); i++) {
            outputFile << ",  " << catalystElementNumbers[i];
        }
    }

    outputFile << "]\n    },\n    \"coords\": {\n      \"3d\": [";

    if(numAtoms() > 0) {
        for(int i = 0; i < layers.size(); i++) {
            if(layers[i].atomPos.size() > 0) {
                for(int j = 0; j < layers[i].atomPos.size(); j++) {
                    outputFile << ((i == 0 && j == 0) ? "  " : "        ");
                    outputFile << layers[i].atomPos[j].x * 10 << ", " << layers[i].atomPos[j].y * 10 << ", " << layers[i].atomPos[j].z * 10;

                    if(i != layers.size() - 1 || j != layers[i].atomPos.size() - 1) {
                        outputFile << ",\n";
                    }
                }
            }
        }

        for(int i = 0; i < catalystAtomPos.size(); i++) {
            outputFile << ",\n        " << catalystAtomPos[i].x * 10 << ", " << catalystAtomPos[i].y * 10 << ", " << catalystAtomPos[i].z * 10;
        }

        outputFile << " ]\n" << "    }\n  },\n  \"bonds\": {\n    \"connections\": {\n      \"index\": [";
    }



    int offset = 0;

    if(hasBonds()) {
        for(int i = 0; i < layers.size(); i++) {
            for(int j = 0; j < layers[i].bonds.size(); j++) {
                outputFile << ((i == 0 && j == 0) ? "  " : "                 ");
                outputFile << layers[i].bonds[j].first + offset << ", " << layers[i].bonds[j].second + offset;

                if(i != layers.size() - 1 || j != layers[i].bonds.size() - 1) {
                    outputFile << ",\n";
                } else {
                    outputFile << " ]\n    },\n    \"order\": [";
                }
            }

            offset += layers[i].atomPos.size();
        }
    }

    if(offset > 0) {
        for(int i = 0; i < layers.size(); i++) {
            for(int j = 0; j < layers[i].bondOrder.size(); j++) {
                outputFile << " " << layers[i].bondOrder[j];
                if(j == layers[i].bondOrder.size() - 1 && i == layers.size() - 1) {
                    outputFile << " ]\n  }\n}";
                } else {
                    outputFile << ",";
                }
            }
        }
    }

}

void MWNTube::addCatalyst(CrystalStructure structure, double radius, double density, double layerCatalystDist, double regularity) {
    double tolerance = min(structure.a, min(structure.b, structure.c));
    for(int i = 0; i < structure.atomElementNumbers.size(); i++) {
        for(int j = i + 1; j < structure.atomElementNumbers.size(); j++) {
            tolerance = min(tolerance, (structure.atomPos(i) - structure.atomPos(j)).magnitude());
        }
    }

    tolerance *= 0.95;

    double layerDiameter = layers[layers.size() - 1].diameter;
    double layerRadius = layerDiameter/2;
    double layerLength = layers[layers.size() - 1].length;
    double adLength = layerLength - (2 * radius);

    double area = layerDiameter * adLength * M_PI;
    int num = max(density * area, 1.1);
    int lNum = max(sqrt(density) * 50 * adLength / (layerDiameter * M_PI), 1.1);
    double rNum = num/double(lNum);
    srand(time(NULL));
    int boundLength = 2 * int(radius/min(structure.a, min(structure.b, structure.c))) + 1;
    int boundOffset = int(radius/min(structure.a, min(structure.b, structure.c)));

    structure.alpha = structure.alpha * M_PI / 180;
    structure.beta = structure.beta * M_PI / 180;
    structure.gamma = structure.gamma * M_PI / 180;

    vector<Vec3d> centers;
    vector<vector<Vec3d>> atoms;
    vector<vector<int>> elements;



    for(int h = 0; h < num; h++) {
        double angle = (floatMod(rand()*0.0000000194163459, 2 * M_PI) / rNum * (1 - regularity)) + (2 * M_PI * (int(floatMod(h * 1.0, rNum)) + int(h/rNum + 0.0000001)*0.5) / int((int(h / rNum) * rNum) - (rNum * (int(h / rNum) + 1))));
        double zPos;
        if(layerLength > 2 * radius) {
            zPos = (floatMod(rand() * 0.000000019748421, adLength)/(max(1, (lNum - 1) * 2)) * (1 - regularity)) + (radius + (lNum == 1 ? adLength/2 : int(h / rNum) * adLength / (lNum - 1)));
            while(zPos < (radius - 0.0000001) || zPos > (layerLength - radius + 0.0000001)) {
                zPos = (floatMod(rand() * 0.000000019748421, adLength)/(max(1, (lNum - 1) * 2)) * (1 - regularity)) + (radius + (lNum == 1 ? adLength/2 : int(h / rNum) * adLength / (lNum - 1)));
            }
        } else {
            zPos = layerLength/2;
        }



        centers.push_back({layerRadius * cos(angle), layerRadius * sin(angle), zPos});
        atoms.push_back({});
        elements.push_back({});

        vector<bool> centerContact(h, false);

        for(int i = 0; i < boundLength; i++) {
            for(int j = 0; j < boundLength; j++) {
                for(int k = 0; k < boundLength; k++) {
                    for(int l = 0; l < structure.atomElementNumbers.size(); l++) {
                        Vec3d intPos = {i - boundOffset + structure.fractionalAtomPos[l].x, j - boundOffset + structure.fractionalAtomPos[l].y, k - boundOffset + structure.fractionalAtomPos[l].z};

                        Vec3d pos = triclinicToEuclidean(intPos, structure.alpha, structure.beta, structure.gamma, structure.a, structure.b, structure.c);
                        pos = pos + centers[h];
                        bool write = true;

                        if((pos - centers[h]).magSquared() > radius * radius) {
                            write = false;
                        }

                        if(Vec2d{pos.x, pos.y}.magSquared() < (layerRadius + layerCatalystDist) * (layerRadius + layerCatalystDist)) {
                            write = false;
                        }

                        if(write) {
                            for(int i = 0; i < h; i++) {
                                if(centerContact[i]) {
                                    for(int j = 0; j < atoms[i].size(); j++) {
                                        if((pos - atoms[i][j]).magSquared() < tolerance * tolerance) {
                                            write = false;
                                        }
                                    }
                                }
                            }

                            if(!write) {
                                continue;
                            }
                            atoms[h].push_back(pos);
                            elements[h].push_back(structure.atomElementNumbers[l]);
                        }
                    }
                }
            }
        }
    }

    for(int i = 0; i < num; i++) {
        for(int j = 0; j < atoms[i].size(); j++) {
            catalystAtomPos.push_back(atoms[i][j]);
            catalystElementNumbers.push_back(elements[i][j]);
        }
    }
}

class Graphene {
public:
    double bondLength;
    int n;
    int m;
    vector<Vec2d> atomPos;
    vector<pair<int, int>> bonds;
public:
    void addHexToAtoms(pair<int, int> hex);
    void generate();
    void output(ofstream& outputFile);
    Graphene(int n, int m, double bondDistance);
};

Graphene::Graphene(int N, int M, double bondDistance) {
    bondLength = bondDistance;
    n = N;
    m = M;
}

void Graphene::addHexToAtoms(pair<int, int> hexagon) {
    vector<Vec2d> offsets =
        {{ 0                     ,  bondLength  },
         {-sqrt(3) * bondLength/2,  bondLength/2},
         {-sqrt(3) * bondLength/2, -bondLength/2},
         { 0                     , -bondLength  },
         { sqrt(3) * bondLength/2, -bondLength/2},
         { sqrt(3) * bondLength/2,  bondLength/2}};

    for(int i = 0; i < offsets.size(); i++) {
        Vec2d tmp = offsets[i] + Vec2d{(hexagon.first + (hexagon.second/2.0)) * sqrt(3) * bondLength, hexagon.second * (-1.5) * bondLength};
        bool add = true;
        for(int j = 0; j < atomPos.size(); j++) {
            if(atomPos[j].x < tmp.x + 0.000001 && atomPos[j].x > tmp.x - 0.000001 && atomPos[j].y < tmp.y + 0.000001 && atomPos[j].y > tmp.y - 0.000001) {
                add = false;
                break;
            }
        }
        if(add) {
            atomPos.push_back(tmp);
        }
    }
}

void Graphene::generate() {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            addHexToAtoms({i, j});
        }
    }

    for(int i = 0; i < atomPos.size(); i++) {
        for(int j = max(0, i - (6 * m)); j < atomPos.size(); j++) {
            bool add = true;
            if((atomPos[i] - atomPos[j]).magSquared() < (bondLength * bondLength * 1.02) && i != j) {
                for(int k = 0; k < bonds.size(); k++) {
                    if(bonds[k] == pair<int, int>(i, j) || bonds[k] == pair<int, int>(j, i)) {
                        add = false;
                        break;
                    }
                }
                if(add) {
                    bonds.push_back({i, j});
                }
            }
        }
    }
}

void Graphene::output(ofstream& outputFile) {
    outputFile << "{\n  \"chemicalJson\": 1, \n  \"atoms\": {\n    \"elements\": {\n      \"number\": [";

    if(atomPos.size() > 0) {
        for(int i = 0; i < atomPos.size() - 1; i++) {
            outputFile << "  6,";
        }
        outputFile << "  6";
    }

    outputFile << "]\n    },\n    \"coords\": {\n      \"3d\": [";

    if(atomPos.size() > 0) {
        for(int j = 0; j < atomPos.size(); j++) {
            outputFile << ((j == 0) ? "  " : "        ");
            outputFile << atomPos[j].x * 10 << ", " << atomPos[j].y * 10 << ", " << double(0);

            if(j != atomPos.size() - 1) {
                outputFile << ",\n";
            } else {
                outputFile << " ]\n" << "    }\n  },\n  \"bonds\": {\n    \"connections\": {\n      \"index\": [";
            }
        }


    }

    int offset = 0;
    int numBonds = 0;

    if(bonds.size() > 0) {
        for(int j = 0; j < bonds.size(); j++) {
            outputFile << ((j == 0) ? "  " : "                 ");
            outputFile << bonds[j].first + offset << ", " << bonds[j].second + offset;

            if(j != bonds.size() - 1) {
                outputFile << ",\n";
            } else {
                outputFile << " ]\n    },\n    \"order\": [";
            }
        }

        offset += atomPos.size();
        numBonds += bonds.size();

    }

    if(offset > 0) {
        for(int i = 0; i < numBonds - 1; i++) {
            outputFile << " 1,";
        }
        outputFile << "1 ]\n  }\n}";
    }

}

double sp2AtomsSheet(double sp2PerTube, double tubesPerSquareMicrometer, double SquareMicrometers) {
    return sp2PerTube * tubesPerSquareMicrometer * SquareMicrometers;
}

int main()
{
    ofstream outputFile;
    outputFile.open("output.cjson");

//    CNTLayer test{4, 3, 3, 0.1421};
    MWNTube example{8, 8.1, 10, 0.1421, 0.34, false};

//    Graphene test{100, 100, 0.1421};

    CrystalStructure palladium {
        90.0,
        90.0,
        90.0,
        0.38898,
        0.38898,
        0.38898,
        {{0, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}},
        {46, 46, 46, 46}
    };



    example.generate(true);
    example.addCatalyst(palladium, 1, 0.2, 0.15, 1);


    example.output(outputFile);
    outputFile.close();
}
