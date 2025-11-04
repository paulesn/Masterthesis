#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Real_timer.h>
#include <map>
#include <iostream>
#include<fstream>
#include <vector>
#include<cmath>
#include<iomanip>
#include<limits>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Segment_2.h>

#undef NDEBUG
#include<cassert>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_intersections_tag                               Itag;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>  Vertex_base;
typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;



// typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Triangulation_data_structure_2<Vb,Face_base> Tds;

typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Point Point;
typedef CDT::Edge  Edge;
typedef CDT::Finite_vertex_handles Finite_vertex_handles;
typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
// typedef Finite_vertex_handles::iterator         Finite_vertex_handles_iterator;
typedef CDT::Face_handle Face_handle;


using namespace std;

class GlobePoint
{
        public:
        double lat, lon;
        GlobePoint(double _lat, double _lon)
        {
                lat=_lat; lon=_lon;
        }
};
void WGS84toGoogleBing(double lat, double lon, double &x, double &y) {
  x = lon * 20037508.34 / 180;
  y = log(tan((90 + lat) * M_PI / 360)) / (M_PI / 180);
  y = y * 20037508.34 / 180;
}
void GoogleBingtoWGS84Mercator (double x, double y, double &lat, double &lon) {
  lon = (x / 20037508.34) * 180;
  lat = (y / 20037508.34) * 180;

  lat = 180/M_PI * (2 * atan(exp(lat * M_PI / 180)) - M_PI / 2);
}

void toMercator(const double latitude, const double longitude, double &x, double &y, int mapWidth=2000, int mapHeight=1000)
{
        // latitude    = 41.145556; // (φ)
        // longitude   = -73.995;   // (λ)

        //mapWidth    = 200;
        //mapHeight   = 100;

        // get x value
        x = (longitude+180)*(mapWidth/360);

        // convert from degrees to radians
        double latRad = latitude*M_PI/180;

        // get y value
        double mercN = log(tan((M_PI/4)+(latRad/2)));
        y     = (mapHeight/2)-(mapWidth*mercN/(2*M_PI));
}

typedef std::pair<Point_2, double> KillBlob;

void refineTriangulation(const CDT &cdt, std::vector<KillBlob> &toInsert, double longestEdge=10000)
{
  long double lengthSum=0;
  long counter=0;
  double maxRadius2=0, maxEdgeLength2=0, minEdgeLength2=std::numeric_limits<double>::max();
  for (CDT::Face_handle fh : cdt.all_face_handles())
  {
	  if (fh->is_in_domain()&&(!cdt.is_infinite(fh)))	// we have an ocean triangle
	  {
		Edge e0=CDT::Edge(fh,0);
		Edge e1=CDT::Edge(fh,1);
		Edge e2=CDT::Edge(fh,2);
		// check for pure ocean triangle
		bool pureOcean=true;
		if (cdt.is_constrained(CDT::Edge(fh,0)))
			pureOcean=false;
		if (cdt.is_constrained(CDT::Edge(fh,1)))
			pureOcean=false;
		if (cdt.is_constrained(CDT::Edge(fh,2)))
			pureOcean=false;
		if (true)	// only refine pure ocean triangles
		{
			Point p1=e0.first->vertex( (e0.second+1)%3 )->point();
			Point p2=e0.first->vertex( (e0.second+2)%3 )->point();
			Point p3=e1.first->vertex( (e1.second+2)%3 )->point();

			K::Circle_2 myCircle(p1,p2,p3);	
			Point center=myCircle.center();
			double sqradius=(myCircle.squared_radius());
			if (sqradius>maxRadius2)
				maxRadius2=sqradius;
			double edge0len2=K::Segment_2(p1,p2).squared_length();
			double edge1len2=K::Segment_2(p2,p3).squared_length();
			double edge2len2=K::Segment_2(p3,p1).squared_length();

			// idea: do not refine triangles where the longest edge is a boundary
			// as this might lead to even skinnier triangles
			double minLen2=edge0len2;
			double maxLen2=edge0len2;
			bool canRefine=true;
			if (cdt.is_constrained(CDT::Edge(fh,0)))
				canRefine=false;

			if (edge1len2>maxLen2)
			{
				if (cdt.is_constrained(CDT::Edge(fh,1)))
					canRefine=false;
				else
					canRefine=true;
				maxLen2=edge1len2;
			}
			minLen2=min(minLen2,edge1len2);
			if (edge2len2>maxLen2)
			{
				if (cdt.is_constrained(CDT::Edge(fh,2)))
					canRefine=false;
				else
					canRefine=true;
				maxLen2=edge2len2;
			}
			minLen2=min(minLen2,edge2len2);



			if (!cdt.is_constrained(CDT::Edge(fh,0)))
			{	
				counter++;
				lengthSum+=sqrt(edge0len2);
				maxEdgeLength2=max(maxEdgeLength2, edge0len2);
				minEdgeLength2=min(minEdgeLength2, edge0len2);
			}
			if (!cdt.is_constrained(CDT::Edge(fh,1)))
			{	
				counter++;
				lengthSum+=sqrt(edge1len2);
				maxEdgeLength2=max(maxEdgeLength2, edge1len2);
				minEdgeLength2=min(minEdgeLength2, edge1len2);
			}
			if (!cdt.is_constrained(CDT::Edge(fh,2)))
			{	
				counter++;
				lengthSum+=sqrt(edge2len2);
				maxEdgeLength2=max(maxEdgeLength2, edge2len2);
				minEdgeLength2=min(minEdgeLength2, edge2len2);
			}

	//		if ((sqradius>2*minEdge)||(sqradius>10000.0*10000.0))
	//		if ((sqradius>2*minEdge)||(sqradius>90000.0*90000.0))
			if ((canRefine)&&((sqradius>1.5*minLen2)||(sqrt(maxLen2)>longestEdge)))
	//		if (false)
			{
				double lat,lon;
				GoogleBingtoWGS84Mercator(center.x(), center.y(), lat, lon);
				if ( (lat<89)&& (lat>-89) && (lon>-180) && (lon<180))
				{
					toInsert.push_back(KillBlob(center,(sqradius)));
				}
			}
		}
	  }
  }
  cout<<"Max radius2 was: "<<sqrt(maxRadius2)<<" max edge length: "<<sqrt(maxEdgeLength2)<<" min edge length: "<<sqrt(minEdgeLength2)<<"\n";
  cout<<"Average length: "<<lengthSum/counter<<"\n";
}

void pruneSteiner(std::vector<KillBlob> &toPrune)
{
	vector<KillBlob> survivors;
	double xmin=std::numeric_limits<double>::max();
	double xmax=std::numeric_limits<double>::lowest();
	double ymin=xmin;
	double ymax=xmax;
	for(int i=0; i<toPrune.size(); i++)
	{
		xmin=min(xmin, toPrune[i].first.x());
		xmax=max(xmax, toPrune[i].first.x());
		ymin=min(ymin, toPrune[i].first.y());
		ymax=max(ymax, toPrune[i].first.y());
	}
	cout<<"We have for X: "<<xmin<<" to "<<xmax<<" and Y:"<< ymin<<" to "<<ymax<<endl;
	cout<<"Initially with "<<toPrune.size()<<" Steiner points"<<endl;
#define GRIDSIZE 1000
	double xDelta=(xmax-xmin)/GRIDSIZE;
	double yDelta=(ymax-ymin)/GRIDSIZE;
	// we have a one-dimensional vector of vectors of KillBlobs (simulating 2-dim vector)
	vector< vector<KillBlob> > myGrid; 
	for(int i=0; i<GRIDSIZE*GRIDSIZE; i++)
	{
		myGrid.push_back(vector<KillBlob>());
	}
	// sort toPrune according to squared radius?
	std::sort(toPrune.begin(), toPrune.end(), [](auto &left, auto &right) { return left.second > right.second; });
	

	for(int i=0; i<toPrune.size(); i++)
	{
		Point_2 curPoint=toPrune[i].first;
		double sqRadius=toPrune[i].second;
		// check if there is already some other point at distance sqrt(sqRadius)
		int xGrid=min(GRIDSIZE-1,int((curPoint.x()-xmin)/xDelta));
		int yGrid=min(GRIDSIZE-1,int((curPoint.y()-ymin)/yDelta));
		int gridRadius=int(sqrt(sqRadius)/min(xDelta,yDelta))+1;
		bool foundTooClose=false;
		// go over all relevant grid cells and check all points in there
		for(int x=max(0, xGrid-gridRadius); x<= min(GRIDSIZE-1, xGrid+gridRadius); x++)
			for(int y=max(0, yGrid-gridRadius); y<= min(GRIDSIZE-1, yGrid+gridRadius); y++)
			{
				int cellNr=GRIDSIZE*y+x;
				for(int j=0; j<myGrid[cellNr].size(); j++)
				{
					double sqDist=K::Segment_2(curPoint, myGrid[cellNr][j].first).squared_length();
					if ((sqDist<sqRadius)&&(sqDist<myGrid[cellNr][j].second))
						foundTooClose=true;
				}
			}
		if (foundTooClose==false)
		{
			myGrid[GRIDSIZE*yGrid+xGrid].push_back(toPrune[i]);
			survivors.push_back(toPrune[i]);
		}
	}
	toPrune=survivors;
	cout<<"After pruning: "<<toPrune.size()<<endl;
}


void markAllFaces(const CDT &cdt)
{
  for(CDT::Face_handle f : cdt.all_face_handles()){
    f->set_in_domain(false);
  }
}
void
discoverInfiniteComponent(const CDT & ct)
{
  //when this function is called, all faces are set NOT "in_domain"
  Face_handle start = ct.infinite_face();
  std::list<Face_handle> queue;
  queue.push_back(start);
  assert(start->is_in_domain()==false);
  int counter=0;
  while(! queue.empty())
  {
    Face_handle fh = queue.front();
    queue.pop_front();
    if (fh->is_in_domain()==false)
    {
	    fh->set_in_domain(true);
	    counter++;
//	    if (counter%100000==0)
//		    std::cout<<counter<<"--"<<queue.size()<<"\n";

	    for(int i = 0; i < 3; i++)
	    {
	      Face_handle fi = fh->neighbor(i);
	      if((!fi->is_in_domain()) && !ct.is_constrained(CDT::Edge(fh,i)))
	      {
		queue.push_back(fi);
		assert(fi->is_in_domain()==false);
	      }
	    }
    }
  }
}



int main(int argc, char*argv[] )
{


	vector< vector<GlobePoint>> myPolys;

	CGAL::Real_timer myTimer;
	myTimer.start();
	cout<<"Usage: "<<argv[0]<<" <file> [minSize [#refinement rounds] [mayEdgeLength>] \n"; 
	cout<<"\tIMPORTANT: assumes mercator input, generates mercator .graph file and lat/lon .gl file\n";

	if (argc==1)
		exit(0);

	bool nostitch=true;

	ifstream inFile(argv[1]);
	int nofPolys;
	inFile>>nofPolys;
	cout<<"We have "<<nofPolys<<" polygons to read"<<endl;

	int minSize=0;
	if (argc>2)
	{
		minSize=atoi(argv[2]);
	}
	cout<<"We will only keep polygons with at least "<<minSize<<" corners"<<endl;

	int refRounds=0;
	if (argc>3)
		refRounds=atof(argv[3]);
	cout<<"We will perfomr "<<refRounds<<" refinement rounds"<<endl;

	double longestEdge=1000000000;
	if (argc>4)
		longestEdge=atof(argv[4]);
	cout<<"Want to make edges shorter than "<<longestEdge<<endl;


	int outputSize=0;
	int totalSize=0;
	double glLatMin=999, glLonMin=999, glLatMax=-999, glLonMax=-999;
	double xMin=std::numeric_limits<double>::max(), yMin=std::numeric_limits<double>::max(), xMax=std::numeric_limits<double>::lowest(), yMax=std::numeric_limits<double>::lowest();
	vector<int> sizeBucket(30,0);
	for(int i=0; i<nofPolys; i++)
	{
		if (i%100000==0)
			cout<<i<<" "<<flush;
		int polySize;
		inFile>>polySize;
		if (polySize==0)
			cout<<i<<": "<<polySize<<" ";
		vector<GlobePoint> curPoly;
		for(int j=0; j<polySize; j++)
		{
			double lat, lon;
			inFile>>lat>>lon;
			curPoly.push_back(GlobePoint(lat,lon));
			totalSize++;
		}
		int bucket=int(log2(curPoly.size()));
		sizeBucket[bucket]++;
		if ((curPoly.size()>minSize))
		{
			outputSize+=curPoly.size();
			myPolys.push_back(curPoly);
			// update max/min values
			for(int j=0; j<polySize; j++)
			{
				double curLat=myPolys.back()[j].lat;
				double curLon=myPolys.back()[j].lon;
				glLatMin=min(curLat, glLatMin);
				glLatMax=max(curLat, glLatMax);
				glLonMin=min(curLon, glLonMin);
				glLonMax=max(curLon, glLonMax);
				double srcX, srcY;
				// toMercator(myPolys.back()[j].lat, myPolys.back()[j].lon, srcX, srcY);
				//WGS84toGoogleBing(myPolys.back()[j].lat, myPolys.back()[j].lon, srcX, srcY);
				//xMin=min(xMin, srcX); xMax=max(xMax,srcX);
				//yMin=min(yMin, srcY); yMax=max(yMax,srcY);
			}
		}
		//cout<<curPoly.size()<<" log2 of size: "<<int(log2(curPoly.size()))<<endl;
	}
	cout<<endl<<"Total Size: "<<totalSize<<endl;
	cout<<endl<<"Output Size: "<<outputSize<<endl;
	cout<<" size of last poly: "<<myPolys.back().size()<<endl;
	cout<<"Will output "<<outputSize<<" in "<<myPolys.size()<<" polygons "<<endl;
	cout<<"Last latitude of last point: "<<myPolys.back().back().lat<<endl;

	cout<<"**** File input at "<<myTimer.time()<<"s\n";

	cout<<"RANGES: Lat: "<<glLatMin<<" to "<<glLatMax<<"\t Lon: "<<glLonMin<<" to "<<glLonMax<<endl;
	cout<<"RANGES: X: "<<xMin<<" to "<<xMax<<"\t Y: "<<yMin<<" to "<<yMax<<endl;


  CDT cdt;
  std::vector<Point> myPoints;
  std::vector<CDT::Vertex_handle> myVHs;
  std::vector<CDT::Vertex_handle> leftSeamVHs, rightSeamVHs;
  int counter=0;


	// DO CGAL triangulation

	for(int i=0;i<myPolys.size(); i++)
	{
		if (i%10000==0)
			cout<<i<<" "<<flush;

		for(int j=0; j<myPolys[i].size(); j++)
		{
			double curLat=myPolys[i][j].lat;
			double curLon=myPolys[i][j].lon;
			if ((curLon==glLonMin)||(curLon==glLonMax))
			{
				cout<<"cut point: "<<curLat<<" "<<curLon<<endl;
			}
			{
				// toMercator(myPolys[i][j].lat, myPolys[i][j].lon, srcX, srcY);
				//WGS84toGoogleBing(myPolys[i][j].lat, myPolys[i][j].lon, srcX, srcY);
				double srcY=myPolys[i][j].lon, srcX=myPolys[i][j].lat;

				myPoints.push_back(Point(srcX, srcY));
				myVHs.push_back(cdt.insert(myPoints.back()));
				myVHs.back()->info()=counter++;
			}
		}
	}
	std::cout<<"**** Point insertion completed at "<<myTimer.time()<<"s\n";

	
	counter=0;
	for(int i=0;i<myPolys.size(); i++)
	{
			// cout<<"Inserting for poly "<<i<<endl;
			for(int j=0; j<myPolys[i].size(); j++)
			{
				// outFile<<(counter+j)<<" "<<(counter+((j+1)%myPolys[i].size()))<<" 4 2\n";
				cdt.insert_constraint(myVHs[counter+j], myVHs[(counter+((j+1)%myPolys[i].size()))]);

			}
			counter+=myPolys[i].size();
	}
	std::cout<<"**** Segment insertion completed at "<<myTimer.time()<<"s\n";

 
  std::cout<<"Number of vertices:"<<cdt.number_of_vertices()<<" Number of faces: "<<cdt.number_of_faces()<<std::endl;
  assert(cdt.is_valid());
  //CGAL::refine_Delaunay_mesh_2(cdt, Criteria());
  //std::cout<<"Number of vertices after conforming:"<<cdt.number_of_vertices()<<std::endl;



  
  // first set all faces to in_domain

  // std::cout<<"**** infinite component detected at "<<myTimer.time()<<"s\n";

  /// now try to do some refinement
  int oldSize=-1;
  for(int i=0; i<refRounds; i++)
  {
	  cout<<"\n\nRefinement round: "<<i<<" at "<<myTimer.time()<<endl;
	  markAllFaces(cdt);
	  discoverInfiniteComponent(cdt);
  	std::vector<KillBlob> toInsert;
  	refineTriangulation(cdt, toInsert, longestEdge);
  	pruneSteiner(toInsert);
	vector<Point> toInsert2;
	for(int ii=0; ii< toInsert.size(); ii++)
	{
		toInsert2.push_back(toInsert[ii].first);
//		if (i==39)
//			cout<<toInsert[ii].first<<" ";
	}
	cdt.insert(toInsert2.begin(), toInsert2.end());
	cout<<"After insertion of "<<toInsert2.size()<<" Steiner points we have "<<cdt.number_of_vertices()<<endl;
        assert(cdt.is_valid());
	if (cdt.number_of_vertices()==oldSize)
		break;
	else
		oldSize=cdt.number_of_vertices();
  }
  //for(int i=0; i<toInsert.size()/10; i++)
  //	  cdt.insert(toInsert[i]);
  std::cout<<"**** after refinement at "<<myTimer.time()<<"s\n";
  std::cout<<"Number of vertices:"<<cdt.number_of_vertices()<<" Number of faces: "<<cdt.number_of_faces()<<std::endl;

  markAllFaces(cdt);
  std::cout<<"After marking"<<std::endl;
  discoverInfiniteComponent(cdt);
  std::cout<<"**** infinite component for last time detected at "<<myTimer.time()<<"s\n";


  int finVertCount=0;

  // give each finite vertex a unique ID starting from 0 (for later indexing with the edges)
  // THIS RESETS the infos provided above (and now also includes all Steiner vertices)
  // if we have a seam vertex, only set the unique ID for the left SEAM vertex
  // in a second step set the unique ID for the right SEAM vertex to its left version (negative)
  int countSeams=0;
for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it){
   // check wheter the vertex is a stitch vertex
     
	double xcoord=it->point().x(), ycoord=it->point().y();
	double lat, lon;
	GoogleBingtoWGS84Mercator(xcoord, ycoord, lat, lon);
	if ((nostitch==false) && ((lon<-179)||(lon>179)))
	{
		if (lon<-179)
		{
			for(int i=0; i<leftSeamVHs.size(); i++)
				if (leftSeamVHs[i]->point()==it->point())
				{
					countSeams++;
					cout<<it->info()<<"L ";
				}
    			it->info()=finVertCount++;
		}
		else if (lon>179)
		{
			bool isSeam=false;
			for(int i=0; i<rightSeamVHs.size(); i++)
				if (rightSeamVHs[i]->point()==it->point())
				{
					isSeam=true;
					countSeams++;
					cout<<it->info()<<"R ";
				}
    			if (isSeam==false)
				it->info()=finVertCount++;
			else
				it->info()=-1;
		}
	}
	else
    		it->info()=finVertCount++;
  }
  cout<<endl<<"finVertCount="<<finVertCount<<endl;
  cout<<"seamVertCount="<<countSeams<<endl;

  // now let the right seam points refer to the left seam points
  for(int i=0; i<rightSeamVHs.size(); i++)
  {
	  assert(rightSeamVHs[i]->info()==-1);
	  rightSeamVHs[i]->info()=-leftSeamVHs[i]->info();
  }




  // now do not output edges that are shared between two non-marked triangles
  // alternatively: only render triangles where in_domain is set to false (oceans)
  // 	always render edge from smaller to larger node ID

  std::map< std::pair<int,int>, bool> edgeThere;
  std::vector< pair<int,int> > edgesToDraw;

  int countOceanFaces=0;
  int countOceanEdges=0;
  for (CDT::Face_handle fh : cdt.all_face_handles())
  {
	  if (fh->is_in_domain()&&(!cdt.is_infinite(fh)))	// we have an ocean triangle
	  {
		  countOceanFaces++;
		  for(int i=0; i<3; i++)
		  {
			Edge e=CDT::Edge(fh,i);
			if ((!cdt.is_infinite(e.first->vertex( (e.second+1)%3 ))) && (!cdt.is_infinite(e.first->vertex( (e.second+2)%3 ))))
			{
				int i1= e.first->vertex( (e.second+1)%3 )->info();
				int i2= e.first->vertex( (e.second+2)%3 )->info();
				// 
				if (i1<0)
					i1=-i1;
				if (i2<0)
					i2=-i2;
				if (i1>i2)
				{
					int tmp=i1;
					i1=i2;
					i2=tmp;
				}
				pair<int, int> curEdge(i1,i2);
				if (i2>finVertCount)
					cout<<"FUCK: "<<i2<<endl;
				if (edgeThere.count(curEdge)==0)
				{
					edgesToDraw.push_back(curEdge);
					edgeThere[curEdge]=true;
					countOceanEdges++;
				}
			}

		  }
	  }
  }
  std::cout<<"The number of ocean triangles is "<<countOceanFaces<<"\n";
  std::cout<<"The number of ocean edges is "<<countOceanEdges<<"\n";

  int count = 0, total_count=0;;
  for (const Edge& e : cdt.finite_edges())
  {
    if (cdt.is_constrained(e))
      ++count;
    total_count++;
  }
  std::cout << "The number of resulting constrained edges is  ";
  std::cout <<  count << std::endl;
  std::cout << "The total number of edges is  ";
  std::cout <<  total_count << std::endl;


	ofstream outFile(string(argv[1])+string(".gl"));
	outFile<<setprecision(20);
	//outFile<<finVertCount<<"\n"<<total_count<<"\n";
	outFile<<finVertCount<<"\n"<<countOceanEdges<<"\n";


	for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it)
	{
		double xcoord=it->point().x(), ycoord=it->point().y();
		double lat, lon;
		GoogleBingtoWGS84Mercator(xcoord, ycoord, lat, lon);

		if (it->info()>=0)
		{
	    		outFile<< lat <<" "<< lon << "\n";
		}
  	}

	std::cout<<"after writing points"<<std::endl;
	total_count=0;

	for(int i=0; i<edgesToDraw.size(); i++)
	{
		pair<int,int> curEdge=edgesToDraw[i];
		int i1=curEdge.first, i2=curEdge.second;
		outFile<<i1<<" "<<i2<<" 1 1\n";
		assert(i2<finVertCount);
	        total_count++;
	}

	std::cout<<"New count: "<<total_count<<"\n";

	outFile.close();

	std::cout<<"Wrote GL file"<<endl;

	// now also write out als set of triangles
	//
	ofstream outFile2(string(argv[1]+string(".graph")));
	outFile2<<setprecision(20);
	outFile2<<finVertCount<<"\n"<<countOceanFaces<<"\n";


	for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it)
	{
		double xcoord=it->point().x(), ycoord=it->point().y();
		double lat, lon;
		GoogleBingtoWGS84Mercator(xcoord, ycoord, lat, lon);
		if (it->info()>=0)
		{
			outFile2<< xcoord << " " << ycoord << "\n";	
		}
  	}

	std::cout<<"after writing points "<<std::endl;


  for (CDT::Face_handle fh : cdt.all_face_handles())
  {
	  if (fh->is_in_domain())	// we have an ocean triangle
	  {
		  if (!cdt.is_infinite(fh))
		  {
		  	countOceanFaces--;
			Edge e0=CDT::Edge(fh,0);
			Edge e1=CDT::Edge(fh,1);
			int i1= e0.first->vertex( (e0.second+1)%3 )->info();
			int i2= e0.first->vertex( (e0.second+2)%3 )->info();
			int i3= e1.first->vertex( (e1.second+2)%3 )->info();
			if (i1<0)
				i1=-i1;
			if (i2<0)
				i2=-i2;
			if (i3<0)
				i3=-i3;
			outFile2<<i1<<" "<<i2<<" "<<i3<<"\n";
			assert(i1>=0);
			assert(i2>=0);
			assert(i3>=0);
			assert(i1<finVertCount);
			assert(i2<finVertCount);
			assert(i3<finVertCount);

		  }
	  }
  }
  	cout<<"Oceanfaces now: "<<countOceanFaces<<endl;


	outFile2.close();
	cout<<"Wrote GRAPH file"<<endl;

  return 0;
}




