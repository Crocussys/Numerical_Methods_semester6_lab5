import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math


h = 0.1


def animated(i, surface, frames, ax, x, y, t):
    t[0] = round(t[0] + tau, 3)
    ax.clear()
    ax.set(xlabel='X', ylabel='Y', zlabel='Z', title=f"t = {t[0]}c")
    ax.set_zlim(20, 110)
    ax.view_init(azim=225)
    surface = ax.plot_surface(x, y, frames[i], cmap='plasma')
    if (t[0] * 1000) % 100 == 0:
        plt.savefig(f"img\\chart2_{t[0]}.svg", bbox_inches='tight')
    return surface,


def chart(frames):
    x = np.arange(0, lxy + h, h)
    y = np.arange(0, lxy + h, h)
    x, y = np.meshgrid(x, y)
    t = [0]
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surface = ax.plot_surface(x, y, frames[0], cmap='plasma')
    ani = animation.FuncAnimation(fig, animated, len(frames) - 1, fargs=(surface, frames, ax, x, y, t),
                                  interval=tau * 10000, repeat=False)
    plt.show()


def f(_t):
    return 4 * (200 - 200 * _t) * (math.sin(200 * _t) + math.cos(300 * _t))


def explicit(u):
    t = 0
    frames = list()
    while t <= t_max:
        u_new = np.zeros((size, size), float)
        for i in range(size):
            for j in range(1, size - 1):
                if i == 0:
                    u_new[i, j] = u[i, j + 1] - 4 * u[i, j] + u[i, j - 1] + 2 * u[i + 1, j]
                elif i == size - 1:
                    u_new[i, j] = u[i, j + 1] - 4 * u[i, j] + u[i, j - 1] + 2 * u[i - 1, j]
                else:
                    u_new[i, j] = u[i, j + 1] - 4 * u[i, j] + u[i, j - 1] + u[i + 1, j] + u[i - 1, j]
        for i in range(size):
            for j in range(1, size - 1):
                u[i, j] += tau * (c * u_new[i, j] / h ** 2 + f(t))
        frames.append(u.copy())
        print(f"t = {t}")
        t += tau
    chart(frames)


def implicit(u):
    sizee = size ** 2
    t = 0
    frames = list()
    while t <= t_max:
        a = np.full((sizee, sizee), 0, float)
        b = np.full(sizee, 0, float)
        lambd = c * tau / h ** 2
        for i in range(sizee):
            a[i, i] = 1 + 2 * lambd
            if i % size > 1:
                a[i, i - 1] = -lambd
            if i % size < size - 2:
                a[i, i + 1] = -lambd
            b[i] = u[i // size, i % size] + f(t) * tau
            if i % size == 1:
                b[i] += u[i // size, 0] * lambd
            if i % size == size - 2:
                b[i] += u[i // size, size - 1] * lambd
        x = np.linalg.inv(a) @ b
        for i in range(size):
            for j in range(1, size - 1):
                u[i, j] = x[i * size + j]
        frames.append(u.copy())
        print(f"t = {t}")
        t += tau
    chart(frames)


if __name__ == "__main__":
    lxy = float(input("Длина и ширина пластины: "))
    c = float(input("Коэффициент теплопроводности: "))
    t_max = float(input("Сколько секунд расчитать: "))
    size = int(lxy / h + 1)
    u_start = np.full((size, size), 30, dtype=float)
    for ind in range(size):
        u_start[ind, 0] = 30
        u_start[ind, size - 1] = 100
    tau = h ** 2 / 10
    explicit(u_start)
    implicit(u_start)
