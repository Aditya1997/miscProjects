import matplotlib.pyplot as plt
import numpy as np

# Assume example values for the inputs
initial_speed = 29  # m/s
launch_angle = 5.5  # degrees
spin_rate = 1  # radians/s
m = 0.27  # kg (Assuming mass of a volleyball)
starting_height = 2.5  # m
lift_coefficient = 0.25  # Adjust as needed
g = 9.81

def volleyball_trajectory_with_magnus(initial_speed, launch_angle, spin_rate, starting_height, lift_coefficient):
    air_density = 1.293  # Air density (kg/m^3)
    ball_radius = 0.11  # Volleyball radius (m)
    drag_coefficient = 0.5  # Drag coefficient for a sphere

    time_step = 0.01
    time_intervals = np.arange(0, 5, time_step)

    x_positions = []
    y_positions = []
    x = -10  # Start 10 meters behind the net
    y = starting_height  # Set the starting height

    # Calculate initial vertical and horizontal velocities from launch angle and speed
    initial_vx = initial_speed * np.cos(np.radians(launch_angle))
    initial_vy = initial_speed * np.sin(np.radians(launch_angle))

    vx = initial_vx
    vy = initial_vy

    for t in time_intervals:
        ball_total_velocity = np.sqrt(vx**2 + vy**2)

        drag_force_magnitude = 0.5 * air_density * ball_total_velocity**2 * drag_coefficient * np.pi * ball_radius**2

        # Calculate the drag force components in x and y directions
        drag_force_x = -drag_force_magnitude * (vx / ball_total_velocity)
        drag_force_y = -drag_force_magnitude * (vy / ball_total_velocity)

        # Calculate magnus force magnitude
        # spin_factor = spin_rate * ball_radius
        # magnus_force_magnitude = 0.5 * lift_coefficient * np.pi * ball_radius**2 * air_density * spin_factor
        # magnus_force_x = magnus_force_magnitude * (vy / ball_total_velocity)
        # magnus_force_y = magnus_force_magnitude * (vx / ball_total_velocity)

        magnus_force_magnitude = 0.5 * air_density * ball_total_velocity**2 * lift_coefficient * np.pi * ball_radius**2 * spin_rate * ball_radius / ball_total_velocity
        magnus_force_x = magnus_force_magnitude * (vy / ball_total_velocity)
        magnus_force_y = -magnus_force_magnitude * (vx / ball_total_velocity)

        ax = (drag_force_x + magnus_force_x) / m
        ay = -g + (drag_force_y + magnus_force_y) / m

        vx += ax * time_step
        vy += ay * time_step

        x += vx * time_step
        y += vy * time_step

        x_positions.append(x)
        y_positions.append(y)

        if y <= 0 and vy < 0:  # Stop the simulation when the ball touches the ground and is moving downward
            break

    return x_positions, y_positions

# Calculate trajectory with no effects (including gravity and launch angle)
time_step = 0.01
time_intervals = np.arange(0, 10, time_step)
x_positions_no_effects = -10 + initial_speed * np.cos(np.radians(launch_angle)) * time_intervals  # Add initial position here
y_positions_no_effects = starting_height + initial_speed * np.sin(np.radians(launch_angle)) * time_intervals - 0.5 * g * time_intervals**2

# Calculate trajectories with different effects
x_positions_magnus, y_positions_magnus = volleyball_trajectory_with_magnus(initial_speed, launch_angle, spin_rate, starting_height, lift_coefficient)
x_positions_drag_only, y_positions_drag_only = volleyball_trajectory_with_magnus(initial_speed, launch_angle, spin_rate, starting_height, 0)  # Set lift_coefficient to 0 here

# Plot the trajectories
plt.figure(figsize=(10, 6))

# Plot trajectory with both drag and Magnus effect
plt.plot(x_positions_magnus, y_positions_magnus, color='red', label='With Magnus Effect')

# Plot trajectory with drag effect only
plt.plot(x_positions_drag_only, y_positions_drag_only, color='blue', label='With Drag Effect Only')

# Plot trajectory with no effects
plt.plot(x_positions_no_effects, y_positions_no_effects, color='green', label='No Effects')

plt.axhline(y=0, color='black', linestyle='--', label='Ground')
plt.plot([0, 0], [0, 2.43], color='orange', linestyle='--', label='Net (2.43 m)')
plt.axvline(x=-9.144, color='gray', linestyle='--', label='Left End of Court')  # Add this line for the left end
plt.axvline(x=9.144, color='gray', linestyle='--', label='Right End of Court')
plt.title("Volleyball Trajectories")
plt.xlabel("Horizontal Distance (m)")
plt.ylabel("Vertical Distance (m)")
plt.legend()
plt.grid()

# Set y-axis limits to accommodate both above and below ground
plt.xlim([-10, 10])
plt.ylim([-0.5, 7])

plt.tight_layout()
plt.show()
