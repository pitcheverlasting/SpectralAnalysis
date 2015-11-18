
vars = ['Tmax', 'Tmin', 'RHmax', 'RHmin', 'u2', 'uz', "sunshine hours", "cloud", "monthly precipitation"
def Penpan(data, constants, solar, alpha, overest, ...)
    "alpha: Please use a numeric value for the alpha (albedo of evaporative surface)"

    Ta = (Tmax + Tmin)/2
    vs_Tmax = 0.6108 * exp(17.27 * Tmax/(Tmax + 237.3))
    vs_Tmin = 0.6108 * exp(17.27 * Tmin/(Tmin + 237.3))
    vas = (vs_Tmax + vs_Tmin)/2
    vabar = (vs_Tmin * RHmax/100 + vs_Tmax * RHmin/100)/2
    P = 101.3 * ((293 - 0.0065 * constants$Elev)/293)^5.26
    delta = 4098 * (0.6108 * exp((17.27 * Ta)/(Ta + 237.3)))/((Ta +
        237.3)^2)
    gamma = 0.00163 * P/constants$lambda
    d_r2 = 1 + 0.033 * cos(2 * pi/365 * J)
    delta2 = 0.409 * sin(2 * pi/365 * J - 1.39)
    w_s = acos(-tan(constants$lat_rad) * tan(delta2))
    N = 24/pi * w_s
    R_a = (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) *
        sin(delta2) + cos(constants$lat_rad) * cos(delta2) *
        sin(w_s))
    R_so = (0.75 + (2 * 10^-5) * constants$Elev) * R_a
    if (solar == "data") {
        R_s = Rs
    }
    else if (solar != "monthly precipitation") {
        R_s = (constants$as + constants$bs * (n/N)) * R_a
    }
    else {
        R_s = (0.85 - 0.047 * Cd) * R_a
    }
    if (is.null(u2)) {
        u2 = uz * 4.87/log(67.8 * constants$z - 5.42)
    }
    else {
        u2 = u2
    }
    R_nl = constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((Tmax +
        273.2)^4 + (Tmin + 273.2)^4)/2 * (1.35 * R_s/R_so -
        0.35)
    P_rad = 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) *
        (constants$lat)^2
    f_dir = -0.11 + 1.31 * R_s/R_a
    R_span = (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) *
        R_s
    R_npan = (1 - constants$alphaA) * R_span - R_nl
    f_pan_u = 1.201 + 1.621 * u2
    Epenpan.Daily = delta/(delta + constants$ap * gamma) * R_npan/constants$lambda +
        constants$ap * gamma/(delta + constants$ap * gamma) *
            f_pan_u * (vas - vabar)
    if (overest == TRUE) {
        Epenpan.Daily = Epenpan.Daily/1.078
    }
    ET.Daily = Epenpan.Daily
    ET.Monthly = aggregate(ET.Daily, as.yearmon(Date.daily,
        "%m/%y"), FUN = sum)
    ET.Annual = aggregate(ET.Daily, floor(as.numeric(as.yearmon(Date.daily,
        "%m/%y"))), FUN = sum)
    ET.MonthlyAve = ET.AnnualAve = NULL
    for (mon in min(as.POSIXlt(Date.daily)$mon):max(as.POSIXlt(Date.daily)$mon)) {
        i = mon - min(as.POSIXlt(Date.daily)$mon) + 1
        ET.MonthlyAve[i] = mean(ET.Daily[as.POSIXlt(Date.daily)$mon ==
            mon])
    }
    for (year in min(as.POSIXlt(Date.daily)$year):max(as.POSIXlt(Date.daily)$year)) {
        i = year - min(as.POSIXlt(Date.daily)$year) + 1
        ET.AnnualAve[i] = mean(ET.Daily[as.POSIXlt(Date.daily)$year ==
            year])
    }
    ET_formulation = "Penpan"
    # ET_type = "Class-A Pan Evaporation"   ("Evaporative surface: ", Surface)

    Surface = paste("user-defined, albedo =", alpha)
    "Solar radiation data have been used directly for calculating evapotranspiration"
    "Sunshine hour data have been used for calculating incoming solar radiation"
    "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
    "Monthly precipitation data have been used for calculating incoming solar radiation"

    # results = list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly,
    # ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve,
    # ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation,
    # ET_type = ET_type, message1 = message1)

    return
}