import tkinter as tk  # интерфейс(окошко)
import sys  # (система)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import time


class Mayatnik:
    def __init__(self, master):
        self.master = master
        master.title("Математический маятник")
        master.geometry("1300x1300")
        master.configure(bg="#e6b4f0")

        self.g = 9.81
        self.l = 2

        self.alpha = tk.DoubleVar(value=75)  # Угол
        self.lambd = tk.DoubleVar(value=0.2)  # Трение
        self.delta_t = tk.DoubleVar(value=0.01)  # Шаг времени
        self.method = tk.StringVar(value="verlet")

        self.last_zero = None
        self.periods = []  # Периоды прошедшие ноль

        self.master.protocol("WM_DELETE_WINDOW", self.close)

        self.create_widgets()

        self.is_running = False
        self.is_paused = False
        self.pause_frame = -1
        self.start_time = None
        self.ani = None

    def create_widgets(self):
        frame = tk.Frame(self.master)
        frame.configure(bg="#e6b4f0")
        frame.pack(fill="y", expand=True)

        padding = {"padx": 10, "pady": 5}
        entry_width = 15
        button_width = 1

        tk.Label(
            frame, text="Начальный угол(градусы):", font=("Helvetica", 10, "bold")
        ).grid(row=0, column=0, sticky="nsew", **padding)
        tk.Entry(frame, textvariable=self.alpha, width=entry_width).grid(
            row=0, column=1, sticky="nsew", **padding
        )

        tk.Label(
            frame, text="Коэфициент угасания(трения):", font=("Helvetica", 10, "bold")
        ).grid(row=1, column=0, sticky="nsew", **padding)
        tk.Entry(frame, textvariable=self.lambd, width=entry_width).grid(
            row=1, column=1, sticky="nsew", **padding
        )

        tk.Label(frame, text="Время шага:", font=("Helvetica", 10, "bold")).grid(
            row=2, column=0, sticky="nsew", **padding
        )
        tk.Entry(frame, textvariable=self.delta_t, width=entry_width).grid(
            row=2, column=1, sticky="nsew", **padding
        )

        tk.Radiobutton(
            frame,
            text="Верлет",
            font=("Helvetica", 10, "bold"),
            variable=self.method,
            value="verlet",
        ).grid(row=3, column=1, columnspan=1, sticky="nsew", **padding)
        tk.Radiobutton(
            frame,
            text="Эйлер",
            font=("Helvetica", 10, "bold"),
            variable=self.method,
            value="euler",
        ).grid(row=3, column=0, columnspan=1, sticky="nsew", **padding)

        tk.Button(
            frame,
            text="Начать",
            font=("Helvetica", 10, "bold"),
            command=self.start_simulation,
            width=button_width,
        ).grid(row=5, column=0, columnspan=1, sticky="nsew", **padding)
        tk.Button(
            frame,
            text="Остановить",
            font=("Helvetica", 10, "bold"),
            command=self.stop_simulation,
            width=button_width,
        ).grid(row=5, column=1, columnspan=1, sticky="nsew", **padding)
        tk.Button(
            frame,
            text="Прервать",
            font=("Helvetica", 10, "bold"),
            command=self.pause_simulation,
            width=button_width,
        ).grid(row=5, column=2, columnspan=2, sticky="nsew", **padding)

        self.period_label_gugens = tk.Label(
            frame, text="Период по Гюгенсу: Вычисляю...", font=("Helvetica", 10, "bold")
        )
        self.period_label_gugens.grid(
            row=1, column=2, columnspan=2, sticky="nsew", **padding
        )

        self.period_label_integral = tk.Label(
            frame, text="Точный период: Вычисляю...", font=("Helvetica", 10, "bold")
        )
        self.period_label_integral.grid(
            row=2, column=2, columnspan=2, sticky="nsew", **padding
        )

        self.period_label_1 = tk.Label(
            frame,
            text="Период через точку равновесия: Вычисляю...",
            font=("Helvetica", 10, "bold"),
        )
        self.period_label_1.grid(
            row=3, column=2, columnspan=2, sticky="nsew", **padding
        )

        self.time_label = tk.Label(
            frame, text="Время работы: 0с", font=("Helvetica", 10, "bold")
        )
        self.time_label.grid(row=0, column=2, columnspan=2, sticky="nsew", **padding)

        # фигура для анимки
        self.angleg, (self.ax2, self.ax) = plt.subplots(2, 1, figsize=(8, 8))
        self.angleg.set_facecolor("#e6b4f0")

        # Конфигурация subplot для анимации маятника
        self.ax.clear()
        self.ax.set_facecolor("#eafab6")
        (self.line,) = self.ax.plot([], [], "o-", lw=0.5)
        self.ax.set_xlim(-1.2 * self.l, 1.2 * self.l)
        self.ax.set_ylim(-1.2 * self.l, 1.2 * self.l)
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.grid(True)

        # Конфигурация subplot для графика угла-время
        self.ax2.clear()
        self.ax2.set_facecolor("#eafab6")
        (self.line2,) = self.ax2.plot([], [], lw=3, color="purple")
        self.ax2.set_xlabel("Время(секунды)")
        self.ax2.set_ylabel("Угол(радианы)")
        self.ax2.grid(True)

        # Создание холста для фигуры
        self.canvas = FigureCanvasTkAgg(self.angleg, master=frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(
            row=6, column=0, columnspan=4, sticky="nsew", **padding
        )

    def close(self):
        self.stop_simulation()
        self.master.destroy()
        sys.exit()

    def pause_simulation(self):
        if self.is_running and not self.is_paused:
            self.is_paused = True
            if self.ani:
                self.ani.event_source.stop()

    def resume_simulation(self):
        if self.is_running and self.is_paused:
            self.is_paused = False
            if self.ani:
                self.ani.event_source.start()

    def start_simulation(self):

        if self.is_running and not self.is_paused:
            return

        if self.is_paused:
            self.resume_simulation()
            return

        # Начинаем новую симуляцию
        self.is_running = True
        self.is_paused = False
        self.start_time = time.time()

        # Инициализируем параметры симуляции
        alpha = np.radians(self.alpha.get())
        self.angle_time_data = []  # угол в момент времени
        delta_t = self.delta_t.get()
        method = self.method.get()

        # Инициализируем переменные
        n = 10000  # кол-во тиков работы
        self.t = np.linspace(0, n * delta_t, n)
        self.angle = np.zeros(n)  # массив углов
        self.speed = np.zeros(n)  # массив скоростей
        self.angle[0] = alpha
        self.speed[0] = 0  # Начальная скорость
        self.initial_total_energy = None

        if method == "verlet":
            self.verlet_integration(n, delta_t)
        elif method == "euler":
            self.euler_integration(n, delta_t)

        # Запускаем анимацию
        self.ani = FuncAnimation(
            self.angleg,
            self.animate,
            frames=n,
            init_func=self.init_animation,
            blit=True,
            interval=delta_t,
        )

        self.canvas.draw()

    def stop_simulation(self):
        if self.ani is not None:
            self.ani.event_source.stop()
            self.is_running = False
            self.is_paused = False
            self.start_time = None

    # ---------------------------------------------------------------------------------------------------------

    def accel_speed(self, i):  # вычислим ускорение маятника в момент времени i
        return (
            -(self.g / self.l) * np.sin(self.angle[i])
            - self.lambd.get() * self.speed[i]
        )

    def euler_integration(self, n, delta_t):
        for i in range(1, n):
            accel = self.accel_speed(i - 1)  # вычисляем ускорение
            self.speed[i] = (
                self.speed[i - 1] + accel * delta_t
            )  # обновляем скорость маятника
            self.angle[i] = (
                self.angle[i - 1] + self.speed[i - 1] * delta_t
            )  # обновляем угол маятника

    def verlet_integration(self, n, delta_t):

        for i in range(n - 1):
            accel = self.accel_speed(i)  # вычисляем ускорение
            self.angle[i + 1] = (
                self.angle[i] + self.speed[i] * delta_t + 0.5 * accel * delta_t**2
            )  # Обновляем угол на следующем шаге
            # Вычисляем ускорение для следующего шага, используя промежуточное положение и скорость
            a_pr = self.accel_speed(i + 1)
            # Обновляем скорость на следующем шаге, используя среднее ускорение между текущим и следующим шагом
            self.speed[i + 1] = self.speed[i] + 0.5 * (accel + a_pr) * delta_t

    def init_animation(self):  # настройки анимации
        self.ax.set_aspect("equal")
        self.ax2.set_aspect("auto")
        self.ax.set_xlim(-1.2 * self.l, 1.2 * self.l)
        self.ax.set_ylim(-1.2 * self.l, 1.2 * self.l)
        self.line.set_data([], [])
        self.ax2.set_ylim(-np.pi, np.pi)
        self.line2.set_data([], [])
        return self.line, self.line2

    def animate(
        self, i
    ):  # обновление позиций маятника на каждом кадре анимации и данных о периоде

        if self.is_paused:
            return self.line, self.line2

        x = [0, self.l * np.sin(self.angle[i])]
        y = [0, -self.l * np.cos(self.angle[i])]

        self.line.set_data(x, y)
        self.line.set_markerfacecolor("purple")

        current_time = time.time() - self.start_time
        self.time_label.config(text=f"Время работы: {current_time:.2f}с")

        Guigens = self.huygens_period()
        self.period_label_gugens.config(text=f"Период по Гюгенсу: {Guigens:.5f} с")
        periodExact = self.exact_period(max(self.angle))
        self.period_label_integral.config(text=f"Точный период: {periodExact:.5f}с")

        # Добавляем текущий угол и время в список angle_time_data и обновляем график колебаний
        self.angle_time_data.append((self.t[i], self.angle[i]))
        times, angles = zip(*self.angle_time_data)
        self.line2.set_data(times, angles)
        self.ax2.set_xlim(0, self.t[i] + 1)

        # Проверяем, произошло ли пересечение через ноль и обновляем метки с периодом(БАРАБАШКА)
        if self.angle[i] > 0 and self.angle[i + 1] < 0:

            t0 = self.t[i] + (self.angle[i] / (self.angle[i] - self.angle[i + 1])) * (
                self.t[i + 1] - self.t[i]
            )
            if self.last_zero is not None:
                period = t0 - self.last_zero
                self.period_label_1.config(
                    text=f"Период через точку равновесия(0): {period:.5f} s"
                )
                self.periods.append(period)
            self.last_zero = t0

        # Останавливаем анимацию, если угол маятника становится очень маленьким в течение двух последовательных итераций
        if (
            abs(self.angle[i]) < abs(self.angle[0]) * self.delta_t.get() * 0.01
            and abs(self.angle[i + 1]) < abs(self.angle[0]) * self.delta_t.get() * 0.01
        ):
            self.ani.event_source.stop()
            self.is_running = False
            self.print_results()

        return self.line, self.line2

    def exact_period(self, fi_max):  # вычисление точного периода
        return 4.0 * np.sqrt(self.l / self.g) * self.cei1(np.sin(fi_max / 2.0))

    def huygens_period(self):  # # вычисление периода по гюгенсу
        return 2.0 * np.pi * np.sqrt(self.l / self.g)

    def cei1(self, k):  # поиск эпилептического интеграла
        # Преобразование массива k в массив t
        t = 1 - np.square(k)
        a = (
            ((0.01451196212 * t + 0.03742563713) * t + 0.03590092383) * t
            + 0.09666344259
        ) * t + 1.38629436112
        b = (
            ((0.00441787012 * t + 0.03328355346) * t + 0.06880248576) * t
            + 0.12498593597
        ) * t + 0.5

        return a - b * np.log(t)


if __name__ == "__main__":
    root = tk.Tk()
    app = Mayatnik(root)
    root.mainloop()
