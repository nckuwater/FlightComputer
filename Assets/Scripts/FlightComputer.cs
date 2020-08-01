using System;
using ModApi;
using ModApi.Common;
using ModApi.Craft;
using ModApi.Flight.Sim;
using ModApi.Flight.UI;
using ModApi.Ui;
using ModApi.Ui.Inspector;
using ModApi.Mods;
using UnityEngine;


namespace Assets.Scripts{
    public class FlightComputer{
        Vector3d craft_r, target_r, craft_v, target_v, craft_unit_r, target_unit_r;
        double craft_v_value, target_v_value,craft_v_r, target_v_r;
        Vector3d craft_h, target_h;
        double craft_h_value, target_h_value;
        double craft_e, target_e, craft_T, target_T, planet_mu, planet_mass;
        double craft_theta, target_theta;
        double craft_E, target_E, craft_Me, target_Me;
        // advance calculation require value, vector
        double craft_init_theta0, target_init_theta0;
        Vector3d craft_p, craft_q, target_p, target_q;
        Vector3d craft_init_r0, craft_init_v0, target_init_r0, target_init_v0;
        double craft_init_x0, craft_init_y0, target_init_x0, target_init_y0;
        double craft_init_dot_x0, craft_init_dot_y0, target_init_dot_x0, target_init_dot_y0;


        // avoid different PCI transfer.
        string craft_parent_name, target_parent_name; 

        // real time coordinate vector
        Vector3d NavNorth, NavEast, NavR;

        // constant
        const double PI2 = 2*Math.PI;

        public FlightComputer(){
            
        }
        private void OnInitializeCalaulation(){
            InitializeCraftData();
        }

        public void InitializeCraftData(){
            // Init planet data
            InitializePlanetData();
            // will update craft, planet.
            var CraftNode = Game.Instance.FlightScene.CraftNode;

            var CraftControl = CraftNode.Controls;


            var CraftScript = CraftNode.CraftScript;
            var CraftData = CraftScript.Data;
            var CraftFlightData = CraftScript.FlightData;

            ICraftOrbitData CraftOrbitData = CraftFlightData.Orbit;
            IOrbit CraftOrbit = CraftNode.Orbit;

            var Craft_input = Game.Instance.Inputs; // IGameInputs

            craft_r = CraftFlightData.Position;
            craft_unit_r = CraftFlightData.PositionNormalized;
            craft_v = CraftFlightData.Velocity;
            craft_v_value = CraftFlightData.VelocityMagnitude;
            craft_v_r = CraftFlightData.VerticalSurfaceVelocity;
            craft_theta = CraftOrbit.TrueAnomaly;
        
            //craft_h = Vector3d.Cross(craft_r, craft_v);// Vector3d
            craft_h = CraftOrbit.AngularMomentum;
            craft_h_value = CraftOrbit.AngularMomentumMag;
            craft_T = CraftOrbitData.Period;
            craft_e = CraftOrbitData.Eccentricity;

            // setup init value, vector for path calculation
            InitializeInitData(craft_r, craft_v, craft_theta, craft_h_value, craft_e, planet_mu, out craft_init_r0, out craft_init_v0, 
                               out craft_init_theta0, out craft_init_x0, out craft_init_y0, out craft_init_dot_x0, out craft_init_dot_y0, 
                               out craft_p, out craft_q);
            Debug.Log($"craft_p {craft_p}");
            Debug.Log($"craft_q {craft_q}");
        }
        public void InitializeTargetData(){
            // Init planet data
            InitializePlanetData();
            // r, v, h, e, p, q
            var CraftNode = Game.Instance.FlightScene.CraftNode;
            var CraftScript = CraftNode.CraftScript;
            var CraftFlightData = CraftScript.FlightData;
            var NavSphereTarget = CraftFlightData.NavSphereTarget;
            if (NavSphereTarget == null){
                Debug.Log("no target selected");
                return;
            }

            IOrbitNode TargetOrbitNode = NavSphereTarget.OrbitNode;
            IOrbit TargetOrbit = TargetOrbitNode.Orbit;

            target_r = NavSphereTarget.Position;
            target_unit_r = Vector3d.Normalize(target_r);
            target_v = NavSphereTarget.Velocity;
            target_v_value = Vector3d.Magnitude(target_v);
            target_v_r = Vector3d.Dot(target_v, target_unit_r);
            target_theta = TargetOrbit.TrueAnomaly;

            target_h = TargetOrbit.AngularMomentum;
            target_h_value = TargetOrbit.AngularMomentumMag;
            target_T = TargetOrbit.Period;
            target_e = TargetOrbit.Eccentricity;

            // Debug.Log(target_r);

            // setup init value, vector for path calculation
            InitializeInitData(target_r, target_v, target_theta, target_h_value, target_e, planet_mu, out target_init_r0, out target_init_v0,
                               out target_init_theta0, out target_init_x0, out target_init_y0, out target_init_dot_x0, out  target_init_dot_y0,
                               out target_p, out target_q);
            Debug.Log("target_p");
            Debug.Log(target_p);
            Debug.Log("target_q");
            Debug.Log(target_q);

        }
        public void InitializeInitData(Vector3d r, Vector3d v, double theta, double h_value, double e, double mu, 
                                       out Vector3d init_r0, out Vector3d init_v0,out double init_theta,
                                       out double init_x0, out double init_y0,
                                       out double init_dot_x0, out double init_dot_y0, out Vector3d p, out Vector3d q){
            init_r0 = r;
            init_v0 = v;
            init_theta = theta;
            init_x0 = Vector3d.Magnitude(init_r0) * Math.Cos(init_theta);
            init_y0 = Vector3d.Magnitude(init_r0) * Math.Sin(init_theta);
            init_dot_x0 = -((mu/h_value)*(Math.Sin(theta)));
            init_dot_y0 = (mu/h_value)*(e + Math.Cos(theta));
            p = ((init_dot_y0/h_value)*init_r0) - (init_y0/h_value)*init_v0;
            q = (init_r0/init_y0) - (init_x0/init_y0)*p;
        }
        public void InitializePlanetData(){
            var PlanetNode = Game.Instance.FlightScene.CraftNode.CraftScript.FlightData.Orbit.Parent;
            var PlanetData = PlanetNode.PlanetData;
            planet_mass = PlanetData.Mass;
            planet_mu = Constants.GravitationConstant * planet_mass;
        }
        public void calculate_data_in_theta(double theta, double init_theta0, Vector3d init_r0, Vector3d init_v0, 
                                            double h_value, double e, double T,
                                            out Vector3d calc_r, out Vector3d calc_v){
            // be sure do init this
            // theta is the variable to be set.
            
            //double E = calculate_E_in_t(craft_t, craft_T, craft_e);
            // get specific value first (reduce code complexity)
            double init_r0_value = Vector3d.Magnitude(init_r0);
            double delta_theta = theta - init_theta0;

            double calc_r_value = get_r(h_value, mu, e, theta);
            double lagr_f = 1 - ( (mu * calc_r )/Math.Pow(h_value, 2) * (1 - Math.Cos(delta_theta)));
            double lagr_g = ((calc_r * init_r0_value) * Math.Sin(delta_theta)) / h_value;

            // Vector3d calc_r = lagr_f * init_r0 + lagr_g * init_v0;
            calc_r = lagr_f * init_r0 + lagr_g * init_v0;

            double lagr_dot_g = 1 - ( ( (mu * init_r0_value) / Math.Pow(h_value, 2) ) * (1 - Math.Cos(delta_theta)) );
            double lagr_dot_f = (1 / lagr_g) * ( (lagr_f * lagr_dot_g) - 1);
            
            // Vector3d calc_v = lagr_dot_f * init_r0 + lagr_dot_g * init_v0;
            calc_v = lagr_dot_f * init_r0 + lagr_dot_g * init_v0;
        }

        public double GetElapsedTimeBetweenTime(double t1, double t2, double T){
            return (t2 + ( Math.Ceiling( (t1 - t2) / T) ) * T ) - t1;
        }
        public double get_r(double h_value, double mu, double e, double theta){
            return (Math.Pow(h_value, 2)/mu) / (1 + e * Math.Cos(theta));
        }
        public double get_T(double mu, double h_value, double e){
            (PI2/Math.Pow(mu, 2))*Math.Pow((h_value/Math.Sqrt(1-Math.Pow(e, 2))), 3);
        }
        public double get_E(double e, double theta){
            return (Math.Atan2((Math.Sqrt(1 - Math.Pow(e, 2) ) * Math.Sin(theta)), (e + Math.Cos(theta))));
        }
        public double get_Me(double e, double theta){
            return (get_E(e, theta) - e*(Math.Sin(get_E(e, theta))));
        }
        public double get_t(double e, double theta, double T){
            return (get_Me(e, theta)/2*Math.PI)*T;
        }
        public double formula_E_zero(double E, double e, double Me){
            return E - e * Math.Sin(E) - Me;
        }
        public double formula_E_derivative(double E, double e){
            return 1 - e * Math.Cos(E);
        }
        public double calculate_Me_in_t(double t, double T){
            return 2*Math.Pi * t / T;
        }
        public double calculate_E_in_t(double t, double T, double e, double init_E = 0){
            // set init_E with theta can boost
            double Me = calculate_Me_in_t(t, T);
            double E = init_E;
            for(int i=0; i<10; ++i){
                E += formula_E_zero(E, e, Me)/formula_E_derivative(E, e);
            }
            return E;
        }
        public double calculate_theta_in_t(double t, double T, double e, double init_theta = 0){
            // init_theta is the variable, not craft_init_theta0
            double init_E = get_E(e, init_theta);
            double E = calculate_E_in_t(t, T, e, init_E);

        }
        public void updateCoordinateVectors(){
            // update North, East, unit Position.
            // NavNorth, NavEast, NavR
            var CraftNode = Game.Instance.FlightScene.CraftNode;
            var CraftScript = CraftNode.CraftScript;
            var CraftFlightData = CraftScript.FlightData;

            NavNorth = CraftFlightData.North;
            NavEast = CraftFlightData.East;
            NavR = Vector3d.Normalize(CraftFlightData.Position);
        }
        public Vector3d XYZ_to_NER(Vector3d XYZ){
            updateCoordinateVectors();
            return new Vector3d(Vector3d.Dot(XYZ, NavNorth), Vector3d.Dot(XYZ, NavEast), Vector3d.Dot(XYZ, NavR));
        }
        public Vector3d NER_to_XYZ(Vector3d NER){
            updateCoordinateVectors();
            return (NavNorth*NER.x + NavEast*NER.y + NavR*NER.z);
        }
        public double rad2deg(double rad){
            return rad*180/Math.PI;
        }
        public double deg2rad(double deg){
            return deg/180*Math.PI;
        }
        public double NER_to_pitch(Vector3d NER){
            return rad2deg(Math.Atan2(NER.z, Vector3d.Magnitude( new Vector3d(NER.x, NER.y, 0))));
        }
        public double NER_to_Heading(Vector3d NER){
            return rad2deg(Math.Atan2(NER.y, NER.x));
        }
    }
}