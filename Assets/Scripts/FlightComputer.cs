using System;
using ModApi;
using ModApi.Common;
using ModApi.Craft;
using ModApi.Flight.Sim;
using ModApi.Flight.Ui;
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
        Vector3d craft_theta, target_theta, craft_p, craft_q, target_p, target_q;
        Vector3d craft_E, target_E, craft_Me, target_Me;

        string craft_parent_name, target_parent_name; // avoid different PCI transfer.

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
            IOrbite CraftOrbit = CraftNode.Orbit;

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
            target_T = TargetOrbitNode.Period;
            target_e = TargetOrbit.Eccentricity;

            Debug.Log(target_r);
            

        }
        public void InitializePlanetData(){
            var PlanetNode = Game.Instance.FlightScene.CraftNode.CraftScript.FlightData.Orbit.Parent;
        }

        public double GetElapsedTimeBetweenTime(double t1, double t2, double T){
            return (t2 + ( Math.Ceiling( (t1 - t2) / T) ) * T ) - t1;
        }

    }
}