namespace Assets.Scripts
{
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using ModApi;
    using ModApi.Common;
    using ModApi.Ui;
    using ModApi.Ui.Inspector;
    using ModApi.Mods;
    using UnityEngine;

    /// <summary>
    /// A singleton object representing this mod that is instantiated and initialize when the mod is loaded.
    /// </summary>
    public class Mod : ModApi.Mods.GameMod
    {
        /// <summary>
        /// Prevents a default instance of the <see cref="Mod"/> class from being created.
        /// </summary>
        //public FlightComputer flightComputer = new FlightComputer();
        public GameObject ww;
        public FlightComputer flightComputer;
        private Mod() : base()
        {
            flightComputer = (new GameObject("ww")).AddComponent<FlightComputer>();
            Debug.Log("fc comp added");
        }
        
        static Mod(){

        }
        protected override void OnModInitialized() {
            Debug.Log("on mod init");
            // flightComputer = new FlightComputer();
            base.OnModInitialized();

            Game.Instance.UserInterface.AddBuildInspectorPanelAction(
                ModApi.Ui.Inspector.InspectorIds.FlightView,
                OnBuildFlightViewInspectorPanel);
            //flightComputer = ww.AddComponent<FlightComputer>();
            
            //flightComputer.OnInitialized();
            //flightComputer = GameObject.AddComponent<FlightComputer>();            
        }
        private void OnBuildFlightViewInspectorPanel(BuildInspectorPanelRequest request){
            var g = new GroupModel("Flight Computer");
            request.Model.AddGroup(g);

            g.Add(new TextButtonModel(
                "show t_r",
                b => flightComputer.InitializeCraftData()));
            g.Add(new TextButtonModel(
                "e2e calc",
                b => flightComputer.calculate_next_e2e_intercept_launch_time()));
            g.Add(new TextButtonModel(
                "wait test",
                b => flightComputer.start_coroutin()));
        }

        /// <summary>
        /// Gets the singleton instance of the mod object.
        /// </summary>
        /// <value>The singleton instance of the mod object.</value>
        public static Mod Instance { get; } = GetModInstance<Mod>();
    }
}
