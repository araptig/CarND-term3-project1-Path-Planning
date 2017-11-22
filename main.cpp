#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
//#include "state.h"

using namespace std;

const double mps2mph 			= 2.23694;		//meters-per-seconds to mile-per-hour

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi()
{
	return(M_PI);
}

double deg2rad(double x)
{
	return(x * pi() / 180);
}

double rad2deg(double x)
{
	return(x * 180 / pi());
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{//hasData
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos)
  {
    return("");
  }
  else if (b1 != string::npos && b2 != string::npos)
  {
    return s.substr(b1, b2 - b1 + 2);
  }
  return("");
}//hasData

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}
	return(closestWaypoint);
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return(closestWaypoint);
}


vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{// (x,y) --> (s,d)
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{// (s,d) -->  (x,y)
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);
	return {x,y};
}

// ------------------------------------------new code -----------------------------

bool car_in_front(vector<vector<double>> sensor_fusion, double s, double d, double mc_speed)
{// returns a flag is there is a car in front
	bool too_close  =false;
	const bool verbose 		 = true;	//for debug

	const double delta_d_max 	= 1.5; 	// if within this displacement ==> in same lane
	const double delta_t 		= 1; // future time horizon
	const double delta_v        = 1;    // if car in front faster do not worry
	const double min_gap 		= 25;   // threshold on longitudinal distance

    double mc_s = s + (delta_t * mc_speed); 	// main car position on time horizon
	unsigned int i = 0;

    while( (i< sensor_fusion.size()) && (too_close == false))
    {//for each other car
		double oc_d_pos = sensor_fusion[i][6]; //other car's displacement
        if (fabs(d-oc_d_pos) < delta_d_max)
        {//same lane
        	// get other car's velocity
        	double v_x   			= sensor_fusion[i][3];
          	double v_y   			= sensor_fusion[i][4];
          	double oc_speed  = sqrt(v_x*v_x + v_y*v_y);	//velocity

          	// gap
          	double oc_s     	    = sensor_fusion[i][5];
          	oc_s  += (delta_t * oc_speed);  // other car position on time horizon
          	double gap = oc_s - mc_s;

          	// three conditions (1) car is close, (2) car is front, (3) car is slower or close velocity
          	too_close = (oc_s-mc_s < min_gap) && (gap>0) && (oc_speed - mc_speed < delta_v);
          	if(verbose && too_close)
          	{//print
          		cout << "gap = " << gap << ", mc/oc speed=" << mc_speed*mps2mph << " / " << oc_speed*mps2mph << " mph" << endl;
          	}//print
        }//same lane
        i++;
    }//for each other car
	return(too_close);
}// returns a flag is there is a car in front



bool check_lane(vector<vector<double>> sensor_fusion, double s,double d, double speed)
{// returns a flag is there is a car in front
	bool lane_available = true;
	const bool verbose 	= false;		//debug
	const double t_int = 1;    //future time
	const double distance_front = 15;
	const double distance_back  = -10;

	for(unsigned int i=0; i< sensor_fusion.size(); i++)
    {//for each car
		double oc_d_pos = sensor_fusion[i][6]; //other car displacement
        if (fabs(d-oc_d_pos)<2)
        {//detected car in that lane
        	double v_x   			= sensor_fusion[i][3];
          	double v_y   			= sensor_fusion[i][4];
          	double oc_current_s     = sensor_fusion[i][5];
          	double oc_speed  = sqrt(v_x*v_x + v_y*v_y);	//velocity

          	double oc_s  = oc_current_s + (t_int * oc_speed);      //other car's position after delta_t seconds
          	double mc_s = s + (t_int*speed);
          	double current_gap =  oc_current_s - s;
          	double future_gap = oc_s - mc_s;

          	// conditions for n=lane NOT being available
          	bool cond1 = (current_gap < distance_front) && (current_gap>distance_back); //other car in the way now
          	bool cond2 = (future_gap < distance_front) && (future_gap>distance_back);   //other car in the way in future
          	if (cond1 || cond2)
          	{//lane not available
          		lane_available = false;
          		if (verbose)
          			cout << "current gap = " << current_gap << ", future gap=" << future_gap << endl;
          	}//lane not available
        }//detected car in same lane
    }//for each car
	return(lane_available);
}// returns a flag is there is a car in front

double get_d(unsigned int lane, const double width)
{
    return( (double(lane) + 0.5)*width);
}


void XY_map2car(double & _ptsx, double & _ptsy, double _ref_x, double _ref_y, double _ref_yaw)
{//in XY-coordinates: map to car coordinates
	double shift_x = _ptsx - _ref_x;
    double shift_y = _ptsy - _ref_y;

    double cos_p   = cos(_ref_yaw);
    double sin_p   = sin(_ref_yaw);

    _ptsx =  shift_x * cos_p  +  shift_y * sin_p;
    _ptsy = -shift_x * sin_p  +  shift_y * cos_p;
}

vector<double> XY_car2map(double _x,double _y, double _ref_x, double _ref_y, double _ref_yaw)
{
	bool const verbose = false;
	vector<double> _xy_map;

	double cos_p   = cos(_ref_yaw);
	double sin_p   = sin(_ref_yaw);
	_xy_map.push_back(_x*cos_p - _y*sin_p + _ref_x);
	_xy_map.push_back(_x*sin_p + _y*cos_p + _ref_y);
    if(verbose)
    {
    	cout << "car (x,y)=(" << _x<< ", " << _y << "), ";
    	cout << "map (x,y)=(" << _xy_map[0] << ", " << _xy_map[1] << ")" << endl;
    }
    return(_xy_map);
}

int main()
{// main

  // constants
  const string map_file_ 	= "../data/highway_map.csv";
  //const double max_s 		= 6945.554;   	//value before wrapping around the track back to 0
  static const double delta_t           = 0.02;         //sampling interval seconds
  static const double lane_width 		= 4;            //in meters

  // initial car state (mc = main car)
  static unsigned int mc_lane 			= 1;			//main car lane initialization
  static double mc_d_pos = get_d(mc_lane,lane_width);	//main car d displacement initialization
  static double mc_desired_speed		= 0;         	//main car speed initialization (m/s)
  static bool mode_lc 					= false;        //0==>KL, 1==>LC
  static bool first_trajectory 			= false;
  static unsigned int freeze_count 		= 0;			//freeze decisions after LC

  // parameters
  static const double max_speed    		= 48/mps2mph; 	//max car speed meters per second
  static const double max_acc           = 8;			//max acceleration meters per second^2
  static const double delta_v  			= max_acc * delta_t; // change in velocity in one sample
  static const unsigned int freeze_set  = 30;

  static const unsigned int N    		= 60;           //default number of trajectory points
  static const int N_sparse             = 3;            //sparse trajectory points after start
  static const double delta_s_sparse    = 25;    		//(m), incremental distance for sparse spline

  //debug
  static const bool verbose             = false;

  // Waypoint map to read from
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // highway_map --> map_waypoints_*
  ifstream in_map_(map_file_.c_str(), ifstream::in);
  string line;
  while (getline(in_map_, line))
  {// while
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }// while

  uWS::Hub h;
  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx,&map_waypoints_dy]
			   (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode)
  {// "42" at the start of the message means there's a websocket message event.
   // The 4 signifies a websocket message
   // The 2 signifies a websocket event

    if (length && (length > 2) && (data[0] == '4') && (data[1] == '2'))
    {//42
      auto s = hasData(data);

      if (s != "")
      {// s != ""
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry")
        {//event == "telemetry"
         // j[1] is the data JSON object

        	// JSON >> variables
        	// Main car's localization data in map coordinates
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];        // degrees
          	double mc_speed = j[1]["speed"];	 // main car - mph
          	mc_speed /= mps2mph;                 // metric units used throughput

          	//state.update_main_car_state(car_x, car_y, mc_speed);
          	if(verbose)
          	{
          		cout << "loc: (x,y)=(" << car_x << "," << car_y << "), (s,d)=(" << car_s << "," << car_d
          			 << "), yaw=" << car_yaw << " v=" << mc_speed << endl;
          	}
          	// Previous path data given to the Planner (N - processed data)
          	vector<double> previous_path_x = j[1]["previous_path_x"];      		//map coordinates
          	vector<double> previous_path_y = j[1]["previous_path_y"];			//map coordinates
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];
          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];


          	// NEW CODE: find trajectory (next_*_vals) in MAP (X,Y) COORDINATES
          	//-----------------------------------------------------------------

         	if(freeze_count > 0)
         	{
         		cout << "count = " << freeze_count << endl;
         		freeze_count --;
         	}
         	else
         	{//no freeze
         		// check if there is a car in front
         		mc_d_pos	= get_d(mc_lane,lane_width);
         		bool too_close = car_in_front(sensor_fusion, car_s, mc_d_pos, mc_speed);

         		if (too_close==true)
         		{// try to change lane
          			bool lane_0_ok = false;
          			bool lane_2_ok = false;
          			unsigned int current_lane = mc_lane;
          			if (current_lane == 1)
          			{// check lanes 0 & 2
          				// check lane 0
          				lane_0_ok = check_lane(sensor_fusion, car_s, 2.0, mc_speed);
          				if(lane_0_ok ==true)
          				{// move to lane 0
          					mode_lc = true;
          				    mc_lane = 0;
          				    mc_d_pos	= get_d(mc_lane,lane_width);
          				    cout << "switch to lane 0: desired d=" << mc_d_pos << endl;
          				}
          				else
          				{//check lane 2
          					lane_2_ok = check_lane(sensor_fusion, car_s, 10.0, mc_speed);
          			        if(lane_2_ok==true)
          			        {// move to lane 0
          			        	mode_lc = true;
          			            mc_lane = 2;
          			            mc_d_pos	= get_d(mc_lane,lane_width);
          			            cout << "switch to lane 2: desired d=" << mc_d_pos << endl;
          			        }
          			        cout << "lane 0:2 OK" <<  lane_0_ok << ":" << lane_2_ok << endl;
          				}
          			}// current lane is 1
          			else
          			{//check lane 1
          				bool lane_1_ok = check_lane(sensor_fusion, car_s, 6.0, mc_speed);
          				if(lane_1_ok==true)
          				{// move to lane 1
          					mode_lc = true;
          				    mc_lane = 1;
          				    mc_d_pos = get_d(mc_lane,lane_width);
          				    cout << "switch to lane 1" << endl;
          				}
          			}//check lane 1

          			// could not find a switch lane
          			if (mode_lc == false)
          				mc_desired_speed -= delta_v;
          			else
          				freeze_count = freeze_set;
         		}// try to change lane
         		else
         		{// no car in front
         			mc_desired_speed += delta_v;
         			if(mc_desired_speed > max_speed)
         				mc_desired_speed = max_speed;
         		}// no car in front
         	}//no freeze


          	// the goal is to determine next_xy_vals
          	//--------------------------------------
            vector<double> next_x_vals;
            vector<double> next_y_vals;
            {//Trajectory formation
            	// (A) use previous (remaining) path point
            	int samples_to_save = previous_path_x.size();

            	for(unsigned int i=0; i<samples_to_save; i++)
            	{
            		next_x_vals.push_back(previous_path_x[i]);
            		next_y_vals.push_back(previous_path_y[i]);
            	}

            	if( (samples_to_save < N))              // || (mode_lc==true))
            	{//need update
            		// (B) spline formation in XY-car coordinates
            		// - step (1): sparsely spaced  way-points and then
            		// - step (2): spline interpolation
            		vector<double> ptsx;	// sparsely spaced way-points (step 1)
            		vector<double> ptsy;    // sparsely spaced way-points (step 1)
            		double ref_x;
            		double ref_y;
            		double ref_yaw;
            		double ref_s;
            		double ref_d;

            		// step 1.a generate first 2 points of ptsx/y (boundary)
            		if (samples_to_save < 2)
            		{// initialization
            			//starting point is car_x/car_y
            			ref_x = car_x;  					// car coordinates on map
            			ref_y = car_y;  					// car coordinates on map
            			ref_s = car_s;
            			ref_d = car_d;
            			ref_yaw = deg2rad(car_yaw); 		// car orientation on map
            			ptsx.push_back(ref_x-cos(ref_yaw)); // previous
            			ptsy.push_back(ref_y-sin(ref_yaw)); // previous
            			ptsx.push_back(ref_x);
            			ptsy.push_back(ref_y);
            			first_trajectory = true;
            		}// initialization
            		else
            		{// normal operation
            			//starting point is the last previous point
            			ref_x = previous_path_x[samples_to_save-1];  	// car coordinates on map
            			ref_y = previous_path_y[samples_to_save-1];  	// car coordinates on map
            			ref_s = end_path_s;
            			ref_d = end_path_d;
            			double ref_x_prev = previous_path_x[samples_to_save-2];
            			double ref_y_prev = previous_path_y[samples_to_save-2];
            			ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);
            			ptsx.push_back(ref_x_prev);
            			ptsx.push_back(ref_x);
            			ptsy.push_back(ref_y_prev);
            			ptsy.push_back(ref_y);
            		}// normal operation

            		// step 1.b generate additional N_sparse points
            		// load additional points, in XY-map coordinates
            		double d_jump = mc_d_pos - ref_d;
            		double _d;
            		double _s;
            		for(unsigned int i = 0; i < N_sparse; i++)
            		{
            			if(mode_lc)
            			{
            				double t_norm = double(i+1)/double(N_sparse);
            				double t_norm_3 = t_norm*t_norm*t_norm;
            				double t_norm_4 = t_norm*t_norm_3;
            				double t_norm_5 = t_norm*t_norm_4;
            				_d = ref_d + d_jump*(10*t_norm_3 - 15*t_norm_4 + 6*t_norm_5);
            			}
            			else
            			{
            				_d = mc_d_pos;
            			}
            			_s = ref_s + double(i+1)*delta_s_sparse;
            			if (mode_lc)
            				cout << "(" << _s << ", " << _d <<")" << endl;
            			vector<double> xy = getXY(_s,_d,map_waypoints_s, map_waypoints_x, map_waypoints_y); //map: (s,d) --> (x,y)
            			ptsx.push_back(xy[0]);
            			ptsy.push_back(xy[1]);
            		}

            		if (verbose)
            		{
            			cout << endl;
            			cout << "sparse trajectory in XY-map coordinates:" << endl;
            			for(unsigned int i = 0; i < ptsx.size(); i++)
            			{
            				cout << "(x,y)=(" << ptsx[i] << ", " << ptsy[i] << ")" << endl;
            			}
            			cout << endl;
            		}

            		// XY map coordinates --> XY car coordinates
            		for(unsigned int i = 0; i < ptsx.size(); i++)
            		{//map coordinates --> car coordinates
            			XY_map2car(ptsx[i],ptsy[i],ref_x,ref_y,ref_yaw);
            		}

            		if (verbose)
            		{
            			cout << "sparse trajectory in xy-car coordinates:" << endl;
            			for(unsigned int i = 0; i < ptsx.size(); i++)
            			{
            				cout << "(" << ptsx[i] << ", " << ptsy[i] << ")" << endl;
            			}
            			cout << endl;
            		}

            		tk::spline s;	//spline interpolation in xy-car coordinates
            		s.set_points(ptsx,ptsy);
            		{//step 2: populate needed trajectory points using spline interpolation
            			double delta_s     = delta_t * mc_desired_speed;       // incremental distance for trajectory
            			double _x = 0.0;
            			double _y;
            			int ub = N-samples_to_save;
            			//if (mode_lc)
            			//{
            			//	ub = int(t_lc/delta_t);
            			//	cout << "upper bound set to " << ub << endl;
            			//}
            			for(unsigned int i = 0; i < ub; i++)
            			{//i
            				if(first_trajectory)
            				{// x = (1/2) a t^2
            					double _t = double(i)*delta_t;
            					_x = 0.5*max_acc*_t*_t;
            				}
            				else
            				{
            					_x += delta_s; //*cos(ref_yaw);
            				}
            				_y = s(_x);

            				// car coordinates to map coordinates
            				vector<double> xy_map = XY_car2map(_x,_y,ref_x,ref_y,ref_yaw);
            				next_x_vals.push_back(xy_map[0]);
            				next_y_vals.push_back(xy_map[1]);
            			}//i
            			if (first_trajectory)
            			{
            				mc_desired_speed = double(N-1)*delta_t*max_acc; // at the end of trajectory
            				cout << "initialized" << endl;
            				first_trajectory = false;
            			}
            		}//update
            		mode_lc = false;
               }//step 2
               //cout << "speed=" << mc_speed << ", s=" << car_s << ", d="<< car_d << ", prev_size=" << prev_size << endl;
           }//KL
 /*
           {//LC
        	   int N_lc      = int(t_lc/delta_t); 				    // new trajectory of 3 seconds (4/0.02 = 200)
        	   //double N_lc2 = N_lc/2.0;								// half way
        	   double delta_s     = delta_t * mc_desired_speed;    // speed --> displacement: incremental distance for trajectory
        	   //double omega       = 2.5*delta_t; 			       // sigmoid parameter

        	   double d_jump      = mc_d_pos - car_d;
        	   for(unsigned int i = 0; i < N_lc; i++)
        	   {
        		   double _s = car_s + double(i+1)*delta_s;
        		   double t_norm = double(i)/double(N_lc);
        		   double t_norm_3 = t_norm*t_norm*t_norm;
        		   double t_norm_4 = t_norm*t_norm_3;
        		   double t_norm_5 = t_norm*t_norm_4;
        		   double _d = car_d + d_jump*(10*t_norm_3 - 15*t_norm_4 + 6*t_norm_5);
        		   //double _d = car_d + d_jump/(1.0 + exp(-omega*(double(i)-N_lc2)));
        		   cout << "(s, d)=(" << _s << ", " << _d << ")" << endl;
        	       vector<double> xy = getXY(_s,_d,map_waypoints_s, map_waypoints_x, map_waypoints_y); //(s,d) --> (x,y)
        	       next_x_vals.push_back(xy[0]);
        	       next_y_vals.push_back(xy[1]);
        	   }
        	   mode_lc = false;
            }//LC
 */

 /*
         	for(unsigned int i = 0; i < N-prev_size; i++)
         	{
         		double next_s = car_s + double(i+prev_size+1)*delta_s;
         		double next_d = mc_d_pos;
         		// (s,d) --> (x,y)
         		vector<double> xy = getXY(next_s,next_d,map_waypoints_s, map_waypoints_x, map_waypoints_y);
         		next_x_vals.push_back(xy[0]);
         		next_y_vals.push_back(xy[1]);
         		cout << "(s,d)=(" << next_s << "," << next_d << "),  (x,y)=" << xy[0] << "," << xy[1] << ")" << endl;
         	}
         	cout << endl;
*/
         	json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;
          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }//event == "telemetry"
      }//s != ""
      else
      {// manual
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }// manual
    }//42
  });//h.onMessage

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t)
  {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req)
  {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length)
  {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
