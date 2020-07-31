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
        Vector3d craft_theta, target_theta;
        Vector3d craft_E, target_E, craft_Me, target_Me;
        // advance calculation require value, vector
        double craft_init_theta, target_init_theta;
        Vector3d craft_p, craft_q, target_p, target_q;
        Vector3d craft_init_r0, craft_init_v0, target_init_r0, target_init_v0;
        Vector3d craft_init_x0, craft_init_y0, target_init_x0, target_init_y0;
        Vector3d craft_init_dot_x0, craft_init_dot_y0, target_init_dot_x0, target_init_dot_y0;


        // avoid different PCI transfer.
        string craft_parent_name, target_parent_name; 

        // real time coordinate vector
        Vector3d NavNorth, NavEast, NavR;

        public FlightComputer(){
            
        }
        private void OnInitializeCalaulation(){
            InitializeCraftData();
        }

        public void InitializeCraftData(){
            // will update craft, planet.
            var CraftNode = Game.Instance.FlightScene.CraftNode;

            var CraftControl = CraftNode.Controls;


            var CraftScript = CraftNode.CraftScript;
            var CraftData = CraftScript.Data;
            var CraftFlightData = CraftScript.FlightData;

            ICraftOrbitData CraftOrbitData = CraftFlightData.Orbit;
            IOrbit CraftOrbit = CraftNode.Orbit;

            var PlanetNode = CraftOrbitData.Parent;
            var PlanetData = PlanetNode.PlanetData;
            planet_mass = PlanetData.Mass;
            planet_mu = Constants.GravitationConstant * planet_mass;

            var Craft_input = Game.Instance.Inputs; // IGameInputs

            craft_r = CraftFlightData.Position;
            craft_unit_r = CraftFlightData.PositionNormalized;
            craft_v = CraftFlightData.Velocity;
            craft_v_value = CraftFlightData.VelocityMagnitude;
            craft_v_r = CraftFlightData.VerticalSurfaceVelocity;
        
            //craft_h = Vector3d.Cross(craft_r, craft_v);// Vector3d
            craft_h = CraftOrbit.AngularMomentum;
            craft_h_value = CraftOrbit.AngularMomentumMag;
            craft_T = CraftOrbitData.Period;
            craft_e = CraftOrbitData.Eccentricity;

            // setup init value, vector for path calculation

        }
        public void InitializeTargetData(){
            // r, v, h, e, p, q
            var CraftNode = Game.Instance.FlightScene.CraftNode;
            var CraftScript = CraftNode.CraftScript;
            var CraftFlightData = CraftScript.FlightData;
            var NavSphereTarget = CraftFlightData.NavSphereTarget;

            IOrbitNode TargetOrbitNode = NavSphereTarget.OrbitNode;
            IOrbit TargetOrbit = TargetOrbitNode.Orbit;

            target_r = NavSphereTarget.Position;
            target_unit_r = Vector3d.Normalize(target_r);
            target_v = NavSphereTarget.Velocity;
            target_v_value = Vector3d.Magnitude(target_v);
            target_v_r = Vector3d.Dot(target_v, target_unit_r);

            target_h = TargetOrbit.AngularMomentum;
            target_h_value = TargetOrbit.AngularMomentumMag;
            target_T = TargetOrbit.Period;
            target_e = TargetOrbit.Eccentricity;

            Debug.Log(target_r);

            // setup init value, vector for path calculation

        }
        public void InitializePlanetData(){
            var PlanetNode = Game.Instance.FlightScene.CraftNode.CraftScript.FlightData.Orbit.Parent;
        }
        public void calculateCraftPath(){

        }

        public double GetElapsedTimeBetweenTime(double t1, double t2, double T){
            return (t2 + ( Math.Ceiling( (t1 - t2) / T) ) * T ) - t1;
        }
        public double get_E(double e, double theta){
            return (Math.Atan2((Math.Sqrt(1 - e^2 ) * Math.Sin(theta))/(e + Math.Cos(theta))));
        }
        public double get_Me(double e, double theta){
            return (get_E(e, theta) - e*(Math.Sin(get_E(e, theta))));
        }
        public double get_t(double e, double theta, double T){
            return (get_Me(e, theta)/2*Math.PI)*T;
        }
        public updateCoordinateVectors(){
            // update North, East, unit Position.
            // NavNorth, NavEast, NavR
            var CraftNode = Game.Instance.FlightScene.CraftNode;
            var CraftScript = CraftNode.CraftScript;
            var CraftFlightData = CraftScript.FlightData;

            NavNorth = CraftFlightData.North;
            NavEast = CraftFlightData.East;
            NavR = Math.Normalize(CraftFlightData.Position);
        }
        public XYZ_to_NER(Vector3d XYZ){
            updateCoordinateVectors();
            return new Vector3d(Vector3d.Dot(XYZ, NavNorth), Vector3d.Dot(XYZ, NavEast), Vector3d.Dot(XYZ, NavR));
        }
        public NER_to_XYZ(Vector3d NER){
            updateCoordinateVectors();
            return (NavNorth*NER.x + NavEast*NER.y + NavR*NER.z);
        }
        public rad2deg(double rad){
            return rad*180/Math.PI;
        }
        public deg2rad(double deg){
            return deg/180*Math.PI;
        }
        public NER_to_pitch(Vector3d NER){
            rad2deg(NER.z, Vector3d.Magnitude( new Vector3d(NER.x, NER.y, 0)));
        }
        public NER_to_Heading(Vector3d NER){
            rad2deg(Math.Atan2(NER.y, NER.x));
        }
    }
}