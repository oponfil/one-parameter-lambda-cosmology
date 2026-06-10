# -*- coding: utf-8 -*-
"""Статистический тест голографического условия T_Hawking = T_mLambda.

Корректное распространение ошибок с полной ковариацией параметров Planck 2018
(линейное приближение + Монте-Карло), вместо упрощённой оценки
sigma(lnR) = (1/12) * dLambda/Lambda при фиксированном omega_m.

Тестируемая величина:
    R = T_Hawking / T_mLambda,
    T_Hawking = hbar*c / (4*pi*k_B*sqrt(R_Lambda*l_P)),
    T_mLambda = T_0 * (omega_Lambda/omega_m)^(1/3),
    R_Lambda  = c / (H_100 * sqrt(omega_Lambda)).

Ключевое тождество: R зависит только от (omega_m, omega_Lambda, T_0);
явная зависимость от h сокращается: R ~ omega_m^(1/3) * omega_Lambda^(-1/12) / T_0.

Ковариация (ln h, ln Omega_m) восстанавливается из опубликованных маржиналов
H_0, Omega_m и хорошо измеренной комбинации Omega_m h^2; перекрёстная проверка
выполняется по Omega_m h^3 (Planck 2018, ур. (12): 0.09633 +/- 0.00029).
"""

import math

import numpy as np

# --- Фундаментальные константы (CODATA 2018) ---
HBAR = 1.054571817e-34   # J s
C = 2.99792458e8         # m / s
G = 6.67430e-11          # m^3 / (kg s^2)
K_B = 1.380649e-23       # J / K
L_P = math.sqrt(HBAR * G / C**3)          # планковская длина, m
KM_S_MPC = 1.0e3 / 3.0856775814913673e22  # 1 км/с/Мпк в 1/с
H100 = 100.0 * KM_S_MPC                   # 100 км/с/Мпк в 1/с

# --- Реликтовая температура (Fixsen 2009 / Planck 2018) ---
T0_MEAN, T0_SIG = 2.7255, 0.0006  # K

# --- Planck 2018, Table 2 (база LambdaCDM, 68% пределы) ---
DATASETS = {
    "TT,TE,EE+lowE": dict(
        H0=67.27, sH0=0.60, Om=0.3166, sOm=0.0084,
        om_h2=0.1432, s_om_h2=0.0013, om_h3=0.09633, s_om_h3=0.00029),
    "TT,TE,EE+lowE+lensing": dict(
        H0=67.36, sH0=0.54, Om=0.3153, sOm=0.0073,
        om_h2=0.1430, s_om_h2=0.0011, om_h3=0.09633, s_om_h3=0.00030),
    "TT,TE,EE+lowE+lensing+BAO": dict(
        H0=67.66, sH0=0.42, Om=0.3111, sOm=0.0056,
        om_h2=0.14240, s_om_h2=0.00087, om_h3=0.09635, s_om_h3=0.00030),
}

N_MC = 4_000_000
RNG = np.random.default_rng(20260610)


def ln_ratio(h, om, t0):
    """ln(T_Hawking / T_mLambda) для плоской Вселенной: Omega_L = 1 - Omega_m."""
    ol = 1.0 - om
    w_l = ol * h**2
    w_m = om * h**2
    r_lambda = C / (H100 * np.sqrt(w_l))
    t_hawking = HBAR * C / (4.0 * np.pi * K_B * np.sqrt(r_lambda * L_P))
    t_ml = t0 * (w_l / w_m) ** (1.0 / 3.0)
    return np.log(t_hawking / t_ml)


def analyze(name, d):
    h, om = d["H0"] / 100.0, d["Om"]
    sx = d["sH0"] / d["H0"]      # sigma(ln h)
    sy = d["sOm"] / d["Om"]      # sigma(ln Omega_m)
    s_lnwm = d["s_om_h2"] / d["om_h2"]

    # Корреляция rho(ln h, ln Omega_m) из условия var(ln om_h2) = var(y + 2x)
    rho = (s_lnwm**2 - sy**2 - 4.0 * sx**2) / (4.0 * sx * sy)
    rho = max(-1.0, min(1.0, rho))

    # Перекрёстная проверка: предсказываем sigma(Omega_m h^3) = sigma(y + 3x)
    s_lnwm3 = math.sqrt(max(sy**2 + 9 * sx**2 + 6 * rho * sx * sy, 0.0))
    s_om_h3_pred = d["om_h3"] * s_lnwm3

    # --- Центральное значение ---
    lnr0 = float(ln_ratio(h, om, T0_MEAN))

    # --- Линейное распространение: d lnR = a dx + b dy - d lnT0 ---
    ol = 1.0 - om
    a = 0.5
    b = 1.0 / 3.0 + (om / ol) / 12.0
    var_lin = (a**2 * sx**2 + b**2 * sy**2 + 2 * a * b * rho * sx * sy
               + (T0_SIG / T0_MEAN) ** 2)
    sig_lin = math.sqrt(var_lin)

    # --- Монте-Карло (точная нелинейная зависимость) ---
    cov = np.array([[sx**2, rho * sx * sy], [rho * sx * sy, sy**2]])
    xy = RNG.multivariate_normal([math.log(h), math.log(om)], cov, size=N_MC)
    t0_s = RNG.normal(T0_MEAN, T0_SIG, size=N_MC)
    lnr = ln_ratio(np.exp(xy[:, 0]), np.exp(xy[:, 1]), t0_s)
    sig_mc = float(np.std(lnr))
    mean_mc = float(np.mean(lnr))

    # --- Для сравнения: метод исходной рукописи ---
    # sigma(ln Lambda) = sigma(ln omega_Lambda) = sigma(-om/ol * y + 2x)
    g = om / ol
    s_lnL = math.sqrt(g**2 * sy**2 + 4 * sx**2 - 4 * g * rho * sx * sy)
    sig_old = s_lnL / 12.0

    nsig = abs(lnr0) / sig_mc
    p_two = math.erfc(nsig / math.sqrt(2.0))

    print(f"=== {name} ===")
    print(f"  h = {h:.4f}, Omega_m = {om:.4f}, rho(ln h, ln Om) = {rho:+.4f}")
    print(f"  Проверка: sigma(Om h^3) предсказано {s_om_h3_pred:.5f}, "
          f"Planck {d['s_om_h3']:.5f}")
    print(f"  Lambda = {3*(H100*h)**2*ol/C**2:.4e} m^-2  "
          f"(sigma(lnLambda) = {s_lnL*100:.2f}%)")
    t_h = math.exp(lnr0) * T0_MEAN * ((ol) / om) ** (1 / 3)  # T_H = R * T_mL
    t_ml = T0_MEAN * (ol / om) ** (1 / 3)
    print(f"  T_Hawking = {t_h:.4f} K,  T_mLambda = {t_ml:.4f} K")
    print(f"  Отклонение lnR = {lnr0*100:+.3f}%")
    print(f"  sigma(lnR): линейно = {sig_lin*100:.3f}%,  MC = {sig_mc*100:.3f}%"
          f"  (среднее MC {mean_mc*100:+.3f}%)")
    print(f"  Метод рукописи (1/12 * dL/L при fix omega_m): {sig_old*100:.3f}%")
    print(f"  Значимость: {nsig:.2f} sigma  (p = {p_two:.3f}, двусторонний)")
    print()
    return lnr0, sig_mc, nsig


def main():
    print(f"l_P = {L_P:.6e} m,  T_0 = {T0_MEAN} +/- {T0_SIG} K,  "
          f"N_MC = {N_MC:.0e}\n")
    for name, d in DATASETS.items():
        analyze(name, d)


if __name__ == "__main__":
    main()
