import numpy as np
import matplotlib.pyplot as plt

# initialization list
p = 2700    # kg/m^3, rock density
C = 1200    # J/(kg·K), heat capacity
k0 = 1.9    # W/(m·K), coefficient of thermal conductivity
B = 2e-3   # K^(-1), constant for k(T)
T0 = 273    # K, the initial temperature in the rod
T1 = 673    # K, the temperature on the left border
tmax = 7 * 24 * 3600    # sec, time
l = 2   # m, rod length

N = 70  # number of space steps
M = 10000   # number of time steps
dx = l / N  # space step
dt = tmax / M   # time step


def check_stability(p, C, k0, dx, dt):
    stability_condition = (p * C * dx ** 2) / (2 * k0)
    return dt <= stability_condition


def update_temperature_const_k(T, T_new, p, C, k0, dx, dt):
    for i in range(1, N):
        T_new[i] = T[i] + dt * k0 / (p * C) * (T[i + 1] - 2 * T[i] + T[i - 1]) / dx ** 2


def apply_boundary_conditions(T, T1):
    T[0] = T1  # first condition
    T[1] = T[0]  # second
    T[N] = T[N - 1]  # added condition to the right border


def k_temp(T, k0, B, T0):
    return k0 / (1 + B * (T - T0))


def update_temperature_var_k(T, T_new, p, C, k0, B, T0, dx, dt):
    for i in range(1, N):
        k_left = 0.5 * (k_temp(T[i - 1], k0, B, T0) + k_temp(T[i], k0, B, T0))
        k_right = 0.5 * (k_temp(T[i + 1], k0, B, T0) + k_temp(T[i], k0, B, T0))
        T_new[i] = T[i] + dt / (p * C * dx ** 2) * ((k_right * (T[i + 1] - T[i])) - (k_left * (T[i] - T[i - 1])))


def solve_heat_equation(T, T_new, p, C, k0, B, T0, dx, dt, k_case='constant'):
    for n in range(M):
        apply_boundary_conditions(T, T1)

        if k_case == 'constant':
            update_temperature_const_k(T, T_new, p, C, k0, dx, dt)
        elif k_case == 'variable':
            update_temperature_var_k(T, T_new, p, C, k0, B, T0, dx, dt)

        T, T_new = T_new, T

    return T


def plot_temperature_distribution(T, l, N, label):
    x = np.linspace(0, l, N+1)
    plt.plot(x, T, label=label)
    plt.xlabel('Distance x (m)')
    plt.ylabel('Temperature T (K)')
    plt.title('Temperature Distribution Over Time')
    plt.legend()
    plt.show()


T = np.zeros(N+1)
T[:] = T0
T_new = np.zeros_like(T)

if not check_stability(p, C, k0, dx, dt):
    print("dt is too big")
else:
    T_final_constant = solve_heat_equation(T, T_new, p, C, k0, B, T0, dx, dt, k_case='constant')
    plot_temperature_distribution(T_final_constant, l, N, label="Constant k")

    T_final_variable = solve_heat_equation(T, T_new, p, C, k0, B, T0, dx, dt, k_case='variable')
    plot_temperature_distribution(T_final_variable, l, N, label="Variable k")

