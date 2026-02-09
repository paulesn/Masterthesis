#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>

using namespace std;

// returns empty vector if holes exist, and otherwise returns polygons describing the connected components of the passed obstacle restricted to [hMin,hMax[
// assumes obstacle to be non-empty in specified range
vector<vector<pair<int, int>>> getPolygons(int w, int h, int hMin, int hMax, vector<pair<int, int>> obstaclePixels, vector<vector<bool>> isObstacle) {
	// cout << "range: " << hMin << " to " << hMax << ", size: " << obstaclePixels.size() << endl;
	int DIR_LEFT = 0;
	int DIR_DOWN = 1;
	int DIR_RIGHT = 2;
	int DIR_UP = 3;

	vector<vector<pair<int, int>>> polygons;

	vector<bool> row2(w, false);
	vector<vector<bool>> gridIsHandled(h, row2); // true = obstacle which was already handled
	// discover all connected components
	for (int k = 0; k < obstaclePixels.size(); k++) {
		if (obstaclePixels[k].first < hMin || obstaclePixels[k].first >= hMax) {
			continue; // not in vertical slice
		}
		if (gridIsHandled[obstaclePixels[k].first][obstaclePixels[k].second]) {
			continue; // already handled
		}
		// get list of obstacle pixels
		vector<pair<int, int>> obstaclePixelsLocal;
		pair<int, int> firstPixel;
		firstPixel.first = obstaclePixels[k].first;
		firstPixel.second = obstaclePixels[k].second;
		obstaclePixelsLocal.push_back(firstPixel);
		gridIsHandled[firstPixel.first][firstPixel.second] = true;
		int exploredLength = 0;
		while (exploredLength < obstaclePixelsLocal.size()) {
			// check 4 neighbors if within range, add to obstacle if not explored and obstacle
			int rowIdx = obstaclePixelsLocal[exploredLength].first;
			int colIdx = obstaclePixelsLocal[exploredLength].second;
			if (rowIdx > hMin) {
				if (isObstacle[rowIdx-1][colIdx] && !gridIsHandled[rowIdx-1][colIdx]) {
					pair<int, int> nextPixel;
					nextPixel.first = rowIdx-1;
					nextPixel.second = colIdx;
					obstaclePixelsLocal.push_back(nextPixel);
					gridIsHandled[rowIdx-1][colIdx] = true;
				}
			}
			if (rowIdx < hMax - 1) {
				if (isObstacle[rowIdx+1][colIdx] && !gridIsHandled[rowIdx+1][colIdx]) {
					pair<int, int> nextPixel;
					nextPixel.first = rowIdx+1;
					nextPixel.second = colIdx;
					obstaclePixelsLocal.push_back(nextPixel);
					gridIsHandled[rowIdx+1][colIdx] = true;
				}
			}
			if (colIdx > 0) {
				if (isObstacle[rowIdx][colIdx-1] && !gridIsHandled[rowIdx][colIdx-1]) {
					pair<int, int> nextPixel;
					nextPixel.first = rowIdx;
					nextPixel.second = colIdx-1;
					obstaclePixelsLocal.push_back(nextPixel);
					gridIsHandled[rowIdx][colIdx-1] = true;
				}
			}
			if (colIdx < w - 1) {
				if (isObstacle[rowIdx][colIdx+1] && !gridIsHandled[rowIdx][colIdx+1]) {
					pair<int, int> nextPixel;
					nextPixel.first = rowIdx;
					nextPixel.second = colIdx+1;
					obstaclePixelsLocal.push_back(nextPixel);
					gridIsHandled[rowIdx][colIdx+1] = true;
				}
			}
			exploredLength++;
		}

		// identify a top left corner pixel and explore boundary of this connected component from there
		vector<bool> fourFalse(4, false); // top, left, bottom, right is boundary? matches traversal direction numbering.
		vector<vector<bool>> row(w, fourFalse);
		vector<vector<vector<bool>>> boundarySides(h, row);	

		for (int h = 0; h < obstaclePixelsLocal.size(); h++) {
			int rowIdxInit = obstaclePixelsLocal[h].first;
			int colIdxInit = obstaclePixelsLocal[h].second;
			if ((rowIdxInit == hMin || !isObstacle[rowIdxInit-1][colIdxInit]) && (colIdxInit == 0 || !isObstacle[rowIdxInit][colIdxInit-1]) && ((rowIdxInit == hMin) || (colIdxInit == 0) || !isObstacle[rowIdxInit-1][colIdxInit-1])) {
				// begin at top left and walk along boundary, initially downwards, until reaching top left of this pixel again
				// we are always at the top left corner of the pixel identified by rowIdx/colIdx, but this pixel is not necessary on the boundary (or even within bounds!)
				int rowIdx = rowIdxInit;
				int colIdx = colIdxInit;
				int currDir = DIR_DOWN;
				vector<int> directionHistory; // includes one entry for each pixel length of the polygon
				
				// each iteration of this while loop moves along the boundary with a step length of 1 pixel
				while (true) {
					directionHistory.push_back(currDir);
					// step in current direction, mark respective boundary of pixel left of what we just traversed as boundary, decide which direction to continue in next step
					if (currDir == DIR_LEFT) {
						boundarySides[rowIdx][colIdx-1][currDir] = true;
						colIdx--;
						if (rowIdx == rowIdxInit && colIdx == colIdxInit) {
							break; // we are back where we started
						}
						
						// which direction to continue: turn left, go straight, or turn right?
						if (colIdx > 0 && isObstacle[rowIdx][colIdx-1]) { // straight is within bounds and stays with or goes into obstacle
							if (rowIdx > hMin && isObstacle[rowIdx-1][colIdx-1]) { // right is within bounds and follows obstacle
								currDir = (currDir + 3) % 4; // turn right
							} else { // go straight
								currDir = currDir + 0; // go straight
							}
						} else { // straight is out of bounds
							currDir = (currDir + 1) % 4; // turn left
						}
					} else if (currDir == DIR_DOWN) {
						boundarySides[rowIdx][colIdx][currDir] = true;
						rowIdx++;

						// which direction to continue: turn left, go straight, or turn right?
						if (rowIdx < hMax && isObstacle[rowIdx][colIdx]) { // straight is within bounds and stays with or goes into obstacle
							if (colIdx > 0 && isObstacle[rowIdx][colIdx-1]) { // right is within bounds and follows obstacle
								currDir = (currDir + 3) % 4; // turn right
							} else { // go straight
								currDir = currDir + 0; // go straight
							}
						} else { // straight is out of bounds
							currDir = (currDir + 1) % 4; // turn left
						}
					} else if (currDir == DIR_RIGHT) {
						boundarySides[rowIdx-1][colIdx][currDir] = true;
						colIdx++;

						// which direction to continue: turn left, go straight, or turn right?
						if (colIdx < w && isObstacle[rowIdx-1][colIdx]) { // straight is within bounds and stays with or goes into obstacle
							if (rowIdx < hMax && isObstacle[rowIdx][colIdx]) { // right is within bounds and follows obstacle
								currDir = (currDir + 3) % 4; // turn right
							} else { // go straight
								currDir = currDir + 0; // go straight
							}
						} else { // straight is out of bounds
							currDir = (currDir + 1) % 4; // turn left
						}
					} else if (currDir == DIR_UP) {
						boundarySides[rowIdx-1][colIdx-1][currDir] = true;
						rowIdx--;

						// which direction to continue: turn left, go straight, or turn right?
						if (rowIdx > hMin && isObstacle[rowIdx-1][colIdx-1]) { // straight is within bounds and stays with or goes into obstacle
							if (colIdx < w && isObstacle[rowIdx-1][colIdx]) { // right is within bounds and follows obstacle
								currDir = (currDir + 3) % 4; // turn right
							} else { // go straight
								currDir = currDir + 0; // go straight
							}
						} else { // straight is out of bounds
							currDir = (currDir + 1) % 4; // turn left
						}
					}
				}
				
				// compose polygon
				pair<int, int> startPos;
				startPos.first = rowIdxInit;
				startPos.second = colIdxInit;
				vector<pair<int, int>> currPoly;
				currPoly.push_back(startPos);
				pair<int, int> secondPos;
				secondPos.first = rowIdxInit+1;
				secondPos.second = colIdxInit;
				currPoly.push_back(secondPos);

				for (int l = 1; l < directionHistory.size(); l++) {
					if (directionHistory[l] == directionHistory[l-1]) {
						// we are following a straight line - adjust previous corner instead of adding more corners along this line
						if (directionHistory[l] == DIR_LEFT) {
							currPoly[currPoly.size()-1].second = currPoly[currPoly.size()-1].second - 1;
						} else if (directionHistory[l] == DIR_DOWN) {
							currPoly[currPoly.size()-1].first = currPoly[currPoly.size()-1].first + 1;
						} else if (directionHistory[l] == DIR_RIGHT) {
							currPoly[currPoly.size()-1].second = currPoly[currPoly.size()-1].second + 1;
						} else if (directionHistory[l] == DIR_UP) {
							currPoly[currPoly.size()-1].first = currPoly[currPoly.size()-1].first - 1;
						}
					} else {
						if (directionHistory[l] == DIR_LEFT) {
							pair<int, int> nextPos;
							nextPos.first = currPoly[currPoly.size()-1].first;
							nextPos.second = currPoly[currPoly.size()-1].second - 1;
							currPoly.push_back(nextPos);
						} else if (directionHistory[l] == DIR_DOWN) {
							pair<int, int> nextPos;
							nextPos.first = currPoly[currPoly.size()-1].first + 1;
							nextPos.second = currPoly[currPoly.size()-1].second;
							currPoly.push_back(nextPos);
						} else if (directionHistory[l] == DIR_RIGHT) {
							pair<int, int> nextPos;
							nextPos.first = currPoly[currPoly.size()-1].first;
							nextPos.second = currPoly[currPoly.size()-1].second + 1;
							currPoly.push_back(nextPos);
						} else if (directionHistory[l] == DIR_UP) {
							pair<int, int> nextPos;
							nextPos.first = currPoly[currPoly.size()-1].first - 1;
							nextPos.second = currPoly[currPoly.size()-1].second;
							currPoly.push_back(nextPos);
						}
					}
				}
				// remove last entry, because that matches the starting position
				currPoly.pop_back();

				polygons.push_back(currPoly);
				
				break; // we do not need to keep looking for a top left boundary pixel anymore, we just found one
			}
		}

		// explore this connected component: are there any other boundaries? if so, we are not free of holes
		for (int h = 0; h < obstaclePixelsLocal.size(); h++) {
			int rowIdx = obstaclePixelsLocal[h].first;
			int colIdx = obstaclePixelsLocal[h].second;
			if((!boundarySides[rowIdx][colIdx][0] && (rowIdx == hMin || !isObstacle[rowIdx - 1][colIdx])) || // top is not marked as boundary but it actually is a boundary? we have found a hole
				(!boundarySides[rowIdx][colIdx][1] && (colIdx == 0 || !isObstacle[rowIdx][colIdx - 1])) || // left is not marked as boundary but it actually is a boundary? we have found a hole
				(!boundarySides[rowIdx][colIdx][2] && (rowIdx == hMax - 1 || !isObstacle[rowIdx + 1][colIdx])) || // bottom is not marked as boundary but it actually is a boundary? we have found a hole
				(!boundarySides[rowIdx][colIdx][3] && (colIdx == w - 1 || !isObstacle[rowIdx][colIdx + 1]))) { // right is not marked as boundary but it actually is a boundary? we have found a hole
				polygons.clear();
				return polygons;
			}
		}
	} // no holes found in any of the connected components that are made up by obstaclePixels? then we are hole-free!
	return polygons;
}

bool charIsObstacle(char q) {
	/*
		. - passable terrain
		G - passable terrain
		@ - out of bounds
		O - out of bounds
		T - trees (unpassable)
		S - swamp (passable from regular terrain)
		W - water (traversable, but not passable from terrain)
	*/
		char obstacles [] = { '@', 'O', 'T' };
	if (q == '@' || q == 'O' || q == 'T') {
		return true;
	} else {
		return false;
	} 
}


int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "usage: mapToCoastlines <mapFile>" << endl;
		exit(0);
	}

	string line = "";
	ifstream inFile(argv[1]);
	std::getline(inFile, line); // "type octile"
	std::getline(inFile, line); // "height <H>"
	int h = stoi(line.substr(7));
	std::getline(inFile, line); // "width <W>"
	int w = stoi(line.substr(6));
	std::getline(inFile, line); // "map"


	// read map: first line we read is line at the top and has rowIdx 0
	vector<bool> row(w, false);
	vector<vector<bool>> isObstacle(h, row); // true = obstacle
	for (int i = 0; i < h; i++) {
		std::getline(inFile, line);
		for (int j = 0; j < w; j++) {
			if (charIsObstacle(line.at(j))) {
				isObstacle[i][j] = true;
			} else {
				isObstacle[i][j] = false;
			}
		}
	}

	cout << "read grid of size " << w << " " << h << ":" << endl;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			cout << isObstacle[i][j];
		}
		cout << endl;
	}
	cout << "now creating polygons of pixel groups which are connected along sides, and not just through corners" << endl;

	// create polygon(s) for each group of touching obstacle pixels
	vector<vector<pair<int, int>>> polygons;

	vector<bool> row2(w, false);
	vector<vector<bool>> gridIsHandled(h, row2); // true = pixel belongs to obstacle which was already handled
	// discover all connected components
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (isObstacle[i][j] && !gridIsHandled[i][j]) {
				// handle obstacle of which this pixel is a part, by dividing it into vertical slices until it has no holes and then producing polygons for each part

				// get list of obstacle pixels
				vector<pair<int, int>> obstaclePixels;
				pair<int, int> firstPixel;
				firstPixel.first = i;
				firstPixel.second = j;
				obstaclePixels.push_back(firstPixel);
				gridIsHandled[i][j] = true;
				int exploredLength = 0;
				while (exploredLength < obstaclePixels.size()) {
					// check 4 neighbors if within range, add to obstacle if not explored and obstacle
					int rowIdx = obstaclePixels[exploredLength].first;
					int colIdx = obstaclePixels[exploredLength].second;
					if (rowIdx > 0) {
						if (isObstacle[rowIdx-1][colIdx] && !gridIsHandled[rowIdx-1][colIdx]) {
							pair<int, int> nextPixel;
							nextPixel.first = rowIdx-1;
							nextPixel.second = colIdx;
							obstaclePixels.push_back(nextPixel);
							gridIsHandled[rowIdx-1][colIdx] = true;
						}
					}
					if (rowIdx < h - 1) {
						if (isObstacle[rowIdx+1][colIdx] && !gridIsHandled[rowIdx+1][colIdx]) {
							pair<int, int> nextPixel;
							nextPixel.first = rowIdx+1;
							nextPixel.second = colIdx;
							obstaclePixels.push_back(nextPixel);
							gridIsHandled[rowIdx+1][colIdx] = true;
						}
					}
					if (colIdx > 0) {
						if (isObstacle[rowIdx][colIdx-1] && !gridIsHandled[rowIdx][colIdx-1]) {
							pair<int, int> nextPixel;
							nextPixel.first = rowIdx;
							nextPixel.second = colIdx-1;
							obstaclePixels.push_back(nextPixel);
							gridIsHandled[rowIdx][colIdx-1] = true;
						}
					}
					if (colIdx < w - 1) {
						if (isObstacle[rowIdx][colIdx+1] && !gridIsHandled[rowIdx][colIdx+1]) {
							pair<int, int> nextPixel;
							nextPixel.first = rowIdx;
							nextPixel.second = colIdx+1;
							obstaclePixels.push_back(nextPixel);
							gridIsHandled[rowIdx][colIdx+1] = true;
						}
					}
					exploredLength++;
				}

				vector<int> horizontalCuts; // list of horizontal cuts, initially 0 - n
				horizontalCuts.push_back(0);
				horizontalCuts.push_back(h);
				int handledCutIndex = 0; // index of latest cut up to which everything has been divided vertically such that divided polygons contain no holes
				
				while (handledCutIndex < horizontalCuts.size() - 1) {
					if (horizontalCuts[handledCutIndex] == horizontalCuts[handledCutIndex+1]) {
						exit(0);
					}
					// check if obstacle is non-empty in current slice
					bool isEmpty = true;
					for (int k = 0; k < obstaclePixels.size(); k++) {
						if (obstaclePixels[k].first >= horizontalCuts[handledCutIndex] && obstaclePixels[k].first < horizontalCuts[handledCutIndex+1]) {
							isEmpty = false;
							break;
						}
					}
					if (isEmpty) {
						handledCutIndex++;					
					} else {
						// generate polygons for current slice (empty unless there are zero holes remaining)
						vector<vector<pair<int, int>>> slicePolygons = getPolygons(w, h, horizontalCuts[handledCutIndex], horizontalCuts[handledCutIndex+1], obstaclePixels, isObstacle);
						if (slicePolygons.size() == 0) {
							// divide current segment in the middle
							int currHeight = horizontalCuts[handledCutIndex+1] - horizontalCuts[handledCutIndex];
							int newCut = horizontalCuts[handledCutIndex] + currHeight / 2;
							horizontalCuts.push_back(h);
							for (int l = horizontalCuts.size() - 2; l > handledCutIndex + 1; l--) {
								horizontalCuts[l] = horizontalCuts[l-1];
							}
							horizontalCuts[handledCutIndex+1] = newCut;
						} else {
							polygons.insert(polygons.end(), slicePolygons.begin(), slicePolygons.end());
							handledCutIndex++;
						}
					}
				}
			}
		}
	}

	// write out coastlines.txt
	cout << "finished, got " << polygons.size() << " polygons" << endl;
	ofstream outFile(argv[1] + string("-coastlines.txt"));
	outFile << polygons.size() << "\n";
	outFile << setprecision(20);
	for (int i = 0; i < polygons.size(); i++) {
		outFile << polygons[i].size() << "\n";
		for (int j = polygons[i].size() - 1; j >= 0; j--) {
			outFile << polygons[i][j].second << " " << polygons[i][j].first << "\n";
		}
	}
	cout << "written coastlines.txt file" << endl;
}