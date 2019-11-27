#include <string>
#include <cstring>
#include <vector>
#include <cstdio>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <omp.h>

template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;


//std::vector<std::vector<Vector3<double>>> triangles;
//std::vector<Vector3<double>> normals;

    // the plane is represented by (x - _p) /dot _normal = 0
    template <typename T>
    class Plane {
    public:
        Plane(Vector3<T> p, Vector3<T> normal) {
            _p = p;
            _normal = normal;
            _normal.normalize();
        }

        Vector3<T>& p() { return _p; }
        Vector3<T>& normal() { return _normal; }
        
        // return if the point is on plane
        // also fill parameter dist as the signed distance from point to plane
        bool onPlane(Vector3<T> point, T& dist) {
            dist = (point - _p).dot(_normal);
            if (std::fabs(dist) < 1e-6) {
                return true;
            } else {
                return false;
            }
        }

    private:
        Vector3<T> _p;
        Vector3<T> _normal;
    };

    template <typename T>
    class Triangle {
    public:
        Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }

        
        // TODO: HW1
        // part 2.1
        // Implement the function to do intersection between triangle and plane p
        // Input: plane p
        // Output: return pairs for intersections with three edges
        // Hint:
        //      - enumerate three edges of the triangle and do intersection individually
        //      - consider the case that no intersection
        //      - consider how to avoid repeated intersection points in returned list
        std::vector<Vector3<T>> IntersectPlane(Plane<T> p) {
			
            std::vector<Vector3<T>> intersections;
            intersections.clear();
			
			int i = 0;
			int ints_ct = 0;
			while (ints_ct < 2 && i<=2) {
				int i1 = i;
				int i2 = (i + 1)%3; 
				i++;
				Vector3<T> p1 = vertices(i1);
				Vector3<T> p2 = vertices(i2);
				if (std::fabs((p2 - p1).dot(p.normal())) < 1e-6) {
					if (std::fabs((p1 - p.p()).dot(p.normal())) < 1e-6) {
						intersections.push_back(p1);
						intersections.push_back(p2);
						ints_ct += 2;
					}
				}
				else {
					T r1 = (p.p() - p1).dot(p.normal()) / ((p2 - p1).dot(p.normal()));
					if (r1 >= 0 && r1 <= 1) {
						intersections.push_back(p1 + r1 * (p2 - p1));
						ints_ct++;
					}
				}
			}

			if (intersections.size() == 2) {
				Vector3<T> p12 = intersections[0] - intersections[1];
				if (std::sqrt(p12[0] * p12[0] + p12[1] * p12[1] + p12[2] * p12[2]) < 1e-6) {
					intersections.clear();
				}
			}
			
            return intersections;
        }

        // TODO: HW2
        // part 1.1
        // Implement the function to do intersection between triangle and a ray
        // Input: a ray, the ray is represented by an origin position and a direction vector
        // Output: return a real number t, the intersection is origin + dir * t, t = -1 means no intersection
        const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const{
			
			T t = -1.0;
			
			Vector3<T> u = _vertices[1]- _vertices[0];  //triangle edge vector u=V1-V0
			Vector3<T> v= _vertices[2] - _vertices[0];   //triangle edge vector v=V2-V0
			T nx = u[1]*v[2]-u[2]*v[1];
			T ny = -(u[0]*v[2]-u[2]*v[0]);
			T nz = u[0]*v[1]-u[1]*v[0];
			T n_mag = sqrt(nx * nx + ny * ny + nz * nz);
			Vector3<T> n(nx/n_mag,ny/n_mag,nz/n_mag); //unit normal vector of the triangle

			//if not parallel to plane triangle lies in
			T n_dir= (n).dot(dir);
			
			if (std::fabs(n_dir) > -(1e-6)) {
				T r = (n).dot(_vertices[0] - origin) / n_dir;

				//if plane intersect with Ray
				if (r >= -(1e-6)) {
					Vector3<T> P = origin + r * dir; //intersection point of ray and plane

					//determine if intersection point P lies inside or on the triangle
					Vector3<T> w = P - _vertices[0]; //vector w=P-V0
					//Barycentric Coordinate 
					//P=V0+s(V1-V0)+t(V2-V0)=V0+s*u+tt*v
					//w=P-V0=s*u+tt*v
					T denominator = ((u).dot(v))*((u).dot(v))-((u).dot(u))*((v).dot(v));
					T s = ((u.dot(v)) * (v.dot(w)) - (v.dot(v)) * (u.dot(w))) / denominator;
					T tt	= ((u.dot(v)) * (w.dot(u)) - (u.dot(u)) * (w.dot(v))) / denominator;
					if (s >= -(1e-6) && tt >= -(1e-6) && s + tt <= 1+(1e-6)) {
						t = r;
					}
				}
			}
			
			return t;
        }

    private:
        Vector3<T> _vertices[3];
    };

//Read in a .stl file and store the triangles and the triangle normals to 
//"triangles", and "normals"
bool ReadSTL(std::string file_name,
    std::vector<std::vector<Vector3<double>>>& triangles, 
    std::vector<Vector3<double>>& normals) {
        
    FILE* fp = std::fopen(file_name.c_str(), "r");

    if (fp == NULL) {
        printf("No STL file found\n");
        return false;
    }

    triangles.clear();
    normals.clear();

    char input[80];
    for (;;) {
        fscanf(fp, "%s", input);
        if (input == std::string("endsolid")) {
            // reach end of file
            break;
        }
        for (;input != std::string("facet");) {
            fscanf(fp, "%s", input);
        }

        std::vector<Vector3<double>> triangle;
        Vector3<double> normal;
        if (std::is_same<double, float>::value) {
            float nx, ny, nz;
            fscanf(fp, "%s %f %f %f\n", input, &nx, &ny, &nz);
            normal[0] = nx; normal[1] = ny; normal[2] = nz;
        }
        else 
            fscanf(fp, "%s %lf %lf %lf\n", input, &normal[0], &normal[1], &normal[2]);

        fscanf(fp, "%s %s", input, input);

        triangle.clear();
        for (int i = 0;i < 3;++i) {
            Vector3<double> p;
            if (std::is_same<double, float>::value) {
                float px, py, pz;
                fscanf(fp, "%s %f %f %f\n", input, &px, &py, &pz);
                p[0] = px; p[1] = py; p[2] = pz;
            }
            else
                fscanf(fp, "%s %lf %lf %lf\n", input, &p[0], &p[1], &p[2]);
            triangle.push_back(p);
        }
        fscanf(fp, "%s %s", input, input);

        triangles.push_back(triangle);
        normals.push_back(normal);
    }

    fclose(fp);
    return true;
}

//normalize the model so that the model is in [-0.5, 0.5]x[-0.5,0.5]x[-0.5, 0.5]
//and zmin should shift to zmin=-0.5
//light source would then have z=zmid=(zmin+zmax)/2
void normalize_model(std::vector<std::vector<Vector3<double>>>& _triangles){
    double xmin=100000;
    double xmax=-100000;
    double ymin=100000;
    double ymax=-100000;
    double zmin=100000;
    double zmax=-100000;
    
    //loop through all triangles
    for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
        Vector3<double> v0 = _triangles[ii][0];//three vertices of a triangle
        Vector3<double> v1 = _triangles[ii][1];
        Vector3<double> v2 = _triangles[ii][2];
        
        if (v0(0) < xmin) { xmin = v0(0); } 
        if (v1(0) < xmin) { xmin = v1(0); } 
        if (v2(0) < xmin) { xmin = v2(0); }
         
        if (v0(0) > xmax) { xmax = v0(0); } 
        if (v1(0) > xmax) { xmax = v1(0); } 
        if (v2(0) > xmax) { xmax = v2(0); }
        
        if (v0(1) < ymin) { ymin = v0(1); } 
        if (v1(1) < ymin) { ymin = v1(1); } 
        if (v2(1) < ymin) { ymin = v2(1); }
         
        if (v0(1) > ymax) { ymax = v0(1); } 
        if (v1(1) > ymax) { ymax = v1(1); } 
        if (v2(1) > ymax) { ymax = v2(1); }
      
        if (v0(2) < zmin) { zmin = v0(2); } 
        if (v1(2) < zmin) { zmin = v1(2); } 
        if (v2(2) < zmin) { zmin = v2(2); }
        
        if (v0(2) > zmax) { zmax = v0(2); } 
        if (v1(2) > zmax) { zmax = v1(2); } 
        if (v2(2) > zmax) { zmax = v2(2); }
    }
    
    double max_range=xmax-xmin;
    if(ymax-ymin>max_range){max_range=ymax-ymin;}
    if(zmax-zmin>max_range){max_range=zmax-zmin;}
    double xmid=0.5*(xmin+xmax);
    double ymid=0.5*(ymin+ymax);
    
    //loop through all triangles
    for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
        _triangles[ii][0](0)=(_triangles[ii][0](0)-xmid)/max_range;
        _triangles[ii][1](0)=(_triangles[ii][1](0)-xmid)/max_range;
        _triangles[ii][2](0)=(_triangles[ii][2](0)-xmid)/max_range;
        
        _triangles[ii][0](1)=(_triangles[ii][0](1)-ymid)/max_range;
        _triangles[ii][1](1)=(_triangles[ii][1](1)-ymid)/max_range;
        _triangles[ii][2](1)=(_triangles[ii][2](1)-ymid)/max_range;
        
        _triangles[ii][0](2)=(_triangles[ii][0](2)-zmin)/max_range-0.5;
        _triangles[ii][1](2)=(_triangles[ii][1](2)-zmin)/max_range-0.5;
        _triangles[ii][2](2)=(_triangles[ii][2](2)-zmin)/max_range-0.5;
    }
    
}

template <typename T>
class Voxelizer {
public:
    Voxelizer(const std::string& stl_file_name, const T dx=0.01) //dx=0.01 corresponds to about 100*100*100 voxels
        : _dx(dx) {
        // Load triangles from the stl file.
        std::vector<Vector3<T>> normals;
        if (!ReadSTL(stl_file_name, _triangles, normals)) {
                std::cout << "ERROR: cannot read " << stl_file_name << std::endl;
                return;
        }
        normalize_model(_triangles);
        // Compute the bounding box of _triangle and save the results into _pmin.
        _pmin = _triangles[0][0];
        Vector3<T> pmax = _triangles[0][0];
        for (const auto& triangle : _triangles)
                for (const auto& v : triangle) {
                        _pmin = _pmin.cwiseMin(v);
                        pmax = pmax.cwiseMax(v);
                }
        for (int i = 0; i < 3; ++i) {
                _pmin[i] -= _dx;
                pmax[i] += _dx;
        }
        // Compute the number of voxels along each direction.
        for (int i = 0; i < 3; ++i)
                _nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;
        // Initialize the voxel array.
        _voxels = std::vector<std::vector<std::vector<bool>>>(_nvoxel.x(),
                std::vector<std::vector<bool>>(_nvoxel.y(),
                        std::vector<bool>(_nvoxel.z(), false)));
    }

    const Vector3<T> pmin() const { return _pmin; }
    const T dx() const { return _dx; }
    const Vector3<T> pmax() const { return _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx; }
    const Vector3<int> voxel_num() const { return _nvoxel; }


    // TODO: HW2
    // part 2.1.
    // Fill the _voxels array with the correct flag.
    void AdvancedVoxelization() {
        const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2]; // number of voxels in each dimension

        //shoot rays with origins on the y-z grid, in direction of x. (i.e. loop through y and z)
        Vector3<T> ray_dir(1.0, 0.0, 0.0); //direction of the ray
        T dx_half = 0.5 * dx();


        std::vector<T> **intersection = new std::vector<T>*[ny];
        for (int i = 0;i < ny;i++) {
                intersection[i] = new std::vector<T>[nz];
        }

        //loop through all triangles
        for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
            Vector3<T> v0 = _triangles[ii][0];//three vertices of a triangle
            Vector3<T> v1 = _triangles[ii][1];
            Vector3<T> v2 = _triangles[ii][2];
            Triangle<T> currentTria(v0, v1, v2);
            T vy_min = v0[1]; 
            if (v1(1) < vy_min) { vy_min = v1(1); } 
            if (v2(1) < vy_min) { vy_min = v2(1); }
            T vy_max=v0[1]; 
            if (v1(1) > vy_max) { vy_max = v1(1); } 
            if (v2(1) > vy_max) { vy_max = v2(1); }
            T vz_min=v0[2]; 
            if (v1(2) < vz_min) { vz_min = v1(2); } 
            if (v2(2) < vz_min) { vz_min = v2(2); }
            T vz_max=v0[2];  
            if (v1(2) > vz_max) { vz_max = v1(2); } 
            if (v2(2) > vz_max) { vz_max = v2(2); }


            int j_min = (int)((vy_min-_pmin[1])/dx());
            int k_min = (int)((vz_min -_pmin[2])/dx());
            int j_max = (int)((vy_max - _pmin[1])/dx());
            int k_max = (int)((vz_max - _pmin[2])/dx());


            //for each triangle, loop through the rays that can potentially intersect with it
            for(int j=j_min;j<=j_max;j++){
                T center_y = _pmin[1] + dx() * j + dx_half;
                for(int k=k_min;k<=k_max;k++){
                    T center_z = _pmin[2] + dx() * k + dx_half;
                    Vector3<T> ray_o(_pmin[0], center_y, center_z); //origin of the ray
                    T t = currentTria.IntersectRay(ray_o, ray_dir);
                    if (t != -1) {
                            //store the intersection point's x coordinate for the corresponding ray
                            //intersection[j*ny+k].push_back(_pmin[0] + t);
                            intersection[j][k].push_back(_pmin[0] + t);
                    }
                }
            }

        }


        //Loop through the rays
        #pragma omp parallel for num_threads(16)
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {

                //If intersection exists
                if (intersection[j][k].size() != 0) {
                    //sort intersection points's x coordinate
                    //sort edges on a certain layer based on x, ascending order
                    //use Insertion Sorting algorithm
                    //itemsSorted: Number of items that have been sorted so far.
                    for (int itemsSorted = 1; itemsSorted < (signed int)intersection[j][k].size(); itemsSorted++) {
                        // Assume that items A[0], A[1], ... A[itemsSorted-1] 
                        // have already been sorted.  Insert A[itemsSorted]
                        // into the sorted part of the list.

                        T temp = intersection[j][k][itemsSorted];  // The edge to be inserted.
                        int loc = itemsSorted - 1;  // Start at end of list.

                        while (loc >= 0 && intersection[j][k][loc] > temp) {
                                intersection[j][k][loc + 1] = intersection[j][k][loc]; // Bump item from A[loc] up to loc+1.
                                loc = loc - 1;       // Go on to next location.
                        }

                        intersection[j][k][loc + 1] = temp; // Put temp in last vacated space.
                    }
                }

                //determine inside/outside for all the grids along this ray
                for (int i = 0; i < nx; ++i) {

                    if (intersection[j][k].size() == 0) {
                        _voxels[i][j][k] = false;
                    }
                    else {
                        T center_x = _pmin[0] + dx() * i + dx_half;
                        bool search = true; int ii = 0; int ct = 0;
                        while (search && ii < (signed int)intersection[j][k].size() / 2) {
                            if (center_x >= intersection[j][k][2 * ii] && center_x <= intersection[j][k][2 * ii + 1]) {
                                _voxels[i][j][k] = true;
                                search = false;
                            }
                            ii++;
                        }
                        if (search) {
                            _voxels[i][j][k] = false;
                        }

                    }
                }

            }
        }

    }

    double PI=3.14159265358979323846;

    void WriteVoxelToMesh(const std::string& stl_file_name) const {
        const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
        std::vector<std::vector<Vector3<int>>> faces;
        std::vector<Vector3<int>> corners({
            Vector3<int>(0, 0, 0),
            Vector3<int>(0, 0, 1),
            Vector3<int>(0, 1, 0),
            Vector3<int>(0, 1, 1),
            Vector3<int>(1, 0, 0),
            Vector3<int>(1, 0, 1),
            Vector3<int>(1, 1, 0),
            Vector3<int>(1, 1, 1)
        });
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                for (int k = 0; k < nz; ++k) {
                    if (!_voxels[i][j][k]) continue;
                    // Check -x direction.
                    Vector3<int> cmin(i, j, k);
                    if (i == 0 || !_voxels[i - 1][j][k]) {
                        faces.push_back({ cmin + corners[0], cmin + corners[1], cmin + corners[3] });
                        faces.push_back({ cmin + corners[0], cmin + corners[3], cmin + corners[2] });
                    }
                    if (i == nx - 1 || !_voxels[i + 1][j][k]) {
                        faces.push_back({ cmin + corners[4], cmin + corners[6], cmin + corners[7] });
                        faces.push_back({ cmin + corners[4], cmin + corners[7], cmin + corners[5] });
                    }
                    if (j == 0 || !_voxels[i][j - 1][k]) {
                        faces.push_back({ cmin + corners[0], cmin + corners[4], cmin + corners[5] });
                        faces.push_back({ cmin + corners[0], cmin + corners[5], cmin + corners[1] });
                    }
                    if (j == ny - 1 || !_voxels[i][j + 1][k]) {
                        faces.push_back({ cmin + corners[2], cmin + corners[3], cmin + corners[7] });
                        faces.push_back({ cmin + corners[2], cmin + corners[7], cmin + corners[6] });
                    }
                    if (k == 0 || !_voxels[i][j][k - 1]) {
                        faces.push_back({ cmin + corners[0], cmin + corners[2], cmin + corners[6] });
                        faces.push_back({ cmin + corners[0], cmin + corners[6], cmin + corners[4] });
                    }
                    if (k == nz - 1 || !_voxels[i][j][k + 1]) {
                        faces.push_back({ cmin + corners[5], cmin + corners[7], cmin + corners[3] });
                        faces.push_back({ cmin + corners[5], cmin + corners[3], cmin + corners[1] });
                    }
                }
        std::ofstream fout(stl_file_name);
        fout << "solid vcg" << std::endl;
        for (const auto& f : faces) {
            std::vector<Vector3<T>> p;
            for (const auto& fi : f) {
                Vector3<T> v = _pmin + fi.cast<T>() * _dx;
                p.push_back(v);
            }
            const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
            fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
            fout << "    outer loop" << std::endl;
            for (const auto& v : p) {
                fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
            }
            fout << "    endloop" << std::endl;
            fout << "  endfacet" << std::endl;
        }
        fout << "endsolid vcg" << std::endl;
    }

//private:
    std::vector<std::vector<Vector3<T>>> _triangles;
    T _dx;  // The size of each voxel.
    Vector3<T> _pmin;    // The min and max corner of the bounding box.
    Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.
    std::vector<std::vector<std::vector<bool>>> _voxels;   // True <-> voxel is occupied.
};

bool rayThruVoxel(double x0, double y0, double z0, double x1, double y1, double z1, 
        double vxmin, double vxmax, double vymin, double vymax, double vzmin, double vzmax){
    
    double xs=x0; double xb=x1;
    if(x0>x1){xs=x1; xb=x0;}
    double ys=y0; double yb=y1;
    if(y0>y1){ys=y1; yb=y0;}
    double zs=z0; double zb=z1;
    if(z0>z1){zs=z1; zb=z0;}
    
    //fix ray x to be vxmin, find the point (vxmin, y, z) on ray, and see if it goes through voxel
    if((vxmin>=xs && vxmin<=xb)){
        double y_temp=(y1-y0)/(x1-x0)*(vxmin-x0)+y0;
        double z_temp=(z1-z0)/(x1-x0)*(vxmin-x0)+z0;
        if(y_temp>vymin && y_temp<vymax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray x to be vxmax
    if((vxmax>=xs && vxmax<=xb)){
        double y_temp=(y1-y0)/(x1-x0)*(vxmax-x0)+y0;
        double z_temp=(z1-z0)/(x1-x0)*(vxmax-x0)+z0;
        if(y_temp>vymin && y_temp<vymax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray y to be vymin
    if((vymin>=ys && vymin<=yb)){
        double x_temp=(x1-x0)/(y1-y0)*(vymin-y0)+x0;
        double z_temp=(z1-z0)/(x1-x0)*(x_temp-x0)+z0;
        if(x_temp>vxmin && x_temp<vxmax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray y to be vymax
    if((vymax>=ys && vymax<=yb)){
        double x_temp=(x1-x0)/(y1-y0)*(vymax-y0)+x0;
        double z_temp=(z1-z0)/(x1-x0)*(x_temp-x0)+z0;
        if(x_temp>vxmin && x_temp<vxmax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray z to be vzmin
    if((vzmin>=zs && vzmin<=zb)){
        double x_temp=(x1-x0)/(z1-z0)*(vzmin-z0)+x0;
        double y_temp=(y1-y0)/(x1-x0)*(x_temp-x0)+y0;
        if(x_temp>vxmin && x_temp<vxmax && y_temp>vymin && y_temp<vymax){
            return true;
        }
    }
    //fix ray z to be vzmax
    if((vzmax>=zs && vzmax<=zb)){
        double x_temp=(x1-x0)/(z1-z0)*(vzmax-z0)+x0;
        double y_temp=(y1-y0)/(x1-x0)*(x_temp-x0)+y0;
        if(x_temp>vxmin && x_temp<vxmax && y_temp>vymin && y_temp<vymax){
            return true;
        }
    }
    
    return false;
}

void read_image(const char *fp, int pixel_m, int pixel_n, int **image_grid){
    int i;
    //pixel_m is pixel number in vertical direction; number of rows
    //pixel_n is pixel number in horizontal direction; number of columns
    *image_grid=new int[pixel_m*pixel_n]; 
    //read in silhouette image, with 1-inside, 2-edge, 0-outside of geometry
    std::ifstream file (fp, std::ifstream::in);
    //std::ifstream file ("C:/Users/jiayi/Desktop/Kay/projects/CR_Voro/luckyCloverBinary720_012.txt", std::ifstream::in);
    for (i=0; i<pixel_m; i++) {
        for (int j=0; j<pixel_n; j++){
            file >> (*image_grid)[pixel_n*i+j];
            if((*image_grid)[pixel_n*i+j]==1){
            }
        }
    }
}

int main() {
    double dx_half=0.5*voxelizer._dx;
    
    //image grids: front, back, right, left, up, down
    int *fgrid; //5:3
    int *bgrid; //5:3
    int *rgrid; //5:3
    int *lgrid; //5:3
    int *ugrid; //5:5
    int *dgrid; //5:5
    
    //1. model
    //read in stl
    //scale the model 
    //voxelize model 
printf("A.\n");
    Voxelizer<double> voxelizer("spherical_shell.stl", 0.001);
printf("B.\n");
    voxelizer.AdvancedVoxelization();
printf("C.\n");
    //2. patterns
    //read in projection images binary pixel patterns in .txt (6 surrounding walls)
    int b_pixel_m=600;
    int b_pixel_n=1000;
    read_image("santa_deer.txt", b_pixel_m, b_pixel_n, &bgrid);
    
    /*
    //debug
    //pixel_m is pixel number in vertical direction; number of rows
    //pixel_n is pixel number in horizontal direction; number of columns
    bgrid=new int[b_pixel_m*b_pixel_n]; 
    //read in silhouette image, with 1-inside, 2-edge, 0-outside of geometry
    std::ifstream file("santa_deer.txt", std::ifstream::in);
    for (int i=0; i<b_pixel_m; i++) {
        for (int j=0; j<b_pixel_n; j++){
            file >> bgrid[b_pixel_n*i+j];
            if(bgrid[b_pixel_n*i+j]==1){
            }
        }
    }
    */
    
printf("D.\n");
    //define surrounding walls domain
    double wall_xmin=-2.5; //bgrid 
    double wall_xmax=2.5;  //fgrid
    double wall_xrange=wall_xmax-wall_xmin;
    double wall_ymin=-2.5; //lgrid
    double wall_ymax=2.5;  //rgrid
    double wall_yrange=wall_ymax-wall_ymin;
    double wall_zmin=-0.5; //dgrid
    double wall_zmax=2.5;  //ugrid
    double wall_zrange=wall_zmax-wall_zmin;
    
    //3. putting together
    //define a light source location
    double lsx=0.0; 
    double lsy=0.0; 
    double lsz=-0.3; //may change to zmid instead
    //voxel grid (lsi, lsj, lsk) that the light source is in.
    int lsi=(int)((lsx-voxelizer._pmin[0])/voxelizer._dx);
    int lsj=(int)((lsy-voxelizer._pmin[1])/voxelizer._dx);
    int isk=(int)((lsz-voxelizer._pmin[2])/voxelizer._dx);
    
    int nx = voxelizer._nvoxel[0], ny = voxelizer._nvoxel[1], nz = voxelizer._nvoxel[2];
printf("E.\n");


    //loop through all pattern grids
    //1. bgrid: 
    double x = wall_xmin;
    //backwall: x<lsx, ii from 0 to lsi
    
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<b_pixel_m; i++) {
        for (int j=0; j<b_pixel_n; j++){
            //if light ray
            //match image grids to surrounding walls coordinates
            if(bgrid[b_pixel_n*i+j]==1){
                double y=wall_ymin+wall_yrange/(b_pixel_n-1)*j;
                double z=wall_zmax-wall_zrange/(b_pixel_m-1)*i;
                
                for(int ii=0; ii<=lsi;ii++){
                    double common=voxelizer._pmin[0]+voxelizer._dx*ii+dx_half;
                    for(int subii=0; subii<5;subii++){
                        double xx=common+0.2*voxelizer._dx*subii;
                        double yy=(y-lsy)/(x-lsx)*(xx-lsx)+lsy;
                        double zz=(z-lsz)/(x-lsx)*(xx-lsx)+lsz;
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);
                        voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                    }
                }
                /*
                double y=wall_ymin+wall_yrange/(b_pixel_n-1)*j;
                double z=wall_zmax-wall_zrange/(b_pixel_m-1)*i;
                
                //loop through all voxels
                //if intersect, void the voxel
                for (int ii = 0; ii < nx; ++ii){
                    for (int jj = 0; jj < ny; ++jj){
                        for (int kk = 0; kk < nz; ++kk) {
                            if (voxelizer._voxels[ii][jj][kk]){
                                
                                double vxmin=voxelizer._pmin[0]+voxelizer._dx*ii;
                                double vxmax=voxelizer._pmin[0]+voxelizer._dx*(ii+1);
                                double vymin=voxelizer._pmin[1]+voxelizer._dx*jj;
                                double vymax=voxelizer._pmin[1]+voxelizer._dx*(jj+1);
                                double vzmin=voxelizer._pmin[2]+voxelizer._dx*kk;
                                double vzmax=voxelizer._pmin[2]+voxelizer._dx*(kk+1);
                                if(rayThruVoxel(lsx, lsy, lsz, x, y, z, vxmin, vxmax, vymin, vymax, vzmin, vzmax)){
                                    voxelizer._voxels[ii][jj][kk]=false;
                                }
                            }
                        }
                    }
                }
                */
            }
        }
    }
printf("F.\n");
    //output the new voxel model
    voxelizer.WriteVoxelToMesh("sphere_model_test1.stl");
printf("G.\n");
    return 0;
}


