# Decentralised Clock Synchronisation

***
Asynchronous decentralised time synchronisation protocol for satellite networks. Given a satellite network input, protocol will aim to bring satellites into agreement onto a common shared time. Simulation can include disruptors. 

Excute simulation from: "initialise_ASYNC"

# INPUTS: 
- Network of interest
- Simulation parameters (duration, asynchronous parameters, number of monte carlo simulations, inclusion of delays...etc)
- Plotting directory 
# OUTPUTS: 
- 1 rate and time convergence plot per batch of monte carlo simulation
- 1 average rate and time convergence plot calculated over all monte carlo simulations

Runs with MATLAB 2023b.

For details on code, see the following figures: 

Simulation initialisation: 

![Satellite_Time_Sync_ASYNC-Initialise_ASYNC drawio](https://github.com/user-attachments/assets/b3c70b27-70c5-49a5-80d2-6fe3bbe265d3)

Main script:

![Satellite_Time_Sync_ASYNC-Script_ASYNC drawio](https://github.com/user-attachments/assets/a63a7c31-b36c-4705-b36e-4c13b2521606)

Update Process:

![Satellite_Time_Sync_ASYNC-Update Process drawio](https://github.com/user-attachments/assets/c7d30261-75c5-4396-9997-ad15d6c26112)

Broadcast Process:

![Satellite_Time_Sync_ASYNC-Broadcast Process drawio](https://github.com/user-attachments/assets/db115e42-0b0b-456c-84c6-443a57cb4d16)


***
