# Equations of State

The following equation of state (EOS) models are supported:
* Ideal gas
* Mie-Gr&uuml;neisen EOS
* SESAME tables

## Ideal Gas Law



## Mie-Gr&uuml;neisen EOS

The Mie-Gr&uuml;neisen model uses a reference curve, usually the principal Hugoniot, and a thermal parameter to describe a pressure-volume-energy EOS surface (e.g., see Chapter 5 in {cite}`Forbes2012`). 

The thermal Gr&uuml;neisen parameter $\gamma$ is given by

$$ 
\gamma(V) = V \left[\frac{dP}{dE}\right]_V, \label{eq1}\tag{1}
$$

where $\gamma$ is often assumed to only depend on volume. Note that parentheses are used to denote dependent variables.

In the IM Tool, the functional form for the Gr&uuml;neisen parameter is currently limited to

$$
\gamma(V)= \gamma_0 \left[V/V_0\right]^q. \label{eq2}\tag{2}
$$

The options for the form for the Gr&uuml;neisen parameter will be expanded as new equation of state options are developed.

A word of caution: these simple equations for the shock Hugoniot and the Gr&uuml;neisen parameter are unlikely to be thermodynamically self-consistent over wide ranges of volume. Beware using the Mie-Gr&uuml;neisen model over wide ranges of pressure and density. Non-physical results will be obtained with equations that are not self-consistent.

The Mie-Gr&uuml;neisen equation of state defines a thermal pressure term, $P(V,E)$, at specific volume $V$ and specific internal energy $E$ referenced to a known state at $V$, which in this case is the principal shock Hugoniot with pressure $P_H(V)$ and specific internal energy $E_H(V)$. Then the thermal pressure term is given by 

$$ 
P(V,E) = P_H(V)+\frac{\gamma(V)}{V} \left[E(V)-E_H(V)\right]. \label{eq3}\tag{3}
$$

A single shock must satisfy the Rankine-Hugoniot (R-H) conservation equations between an initial state, $s_0$, and the shock state, $s_H$. The initial state $s_0$ is defined by the state variables: specific volume $V_0$, presure $P_0$, and specific energy $E_0$ (which then specifies the initial specific entropy $S_0$ and temperature $T_0$). Note the tabulated initial value for the Gr&uuml;neisen parameter $\gamma_0$ must correspond to state $s_0$. The R-H energy conservation equation is 

$$ 
E_H(V)=E_0+\left[\frac{P_0+P_H(V)}{2}\right]\left[V_0-V\right]. \label{eq4}\tag{4}
$$

Combining these two equations defines the Mie-Gr&uuml;neisen equation of state surface, which is limited to the volume range where $P_H(V)$ is defined,

$$
P(V,E) & =  P_H(V)+\frac{\gamma(V)}{V}\left[E(V)-E_0\right] - \left[P_0+P_H(V)\right]\frac{\gamma(V)}{2V}[V_0-V]  \\
       & =  \frac{\gamma(V)}{V}\left[E(V)-E_0\right] + P_H(V)\left\{ 1- \frac{\gamma(V)}{2V}\left[V_0-V\right] \right\}
        - \frac{P_0\gamma(V)}{2V}\left[V_0-V\right],\label{eq5} (5)
$$

and the last term is dropped when $P_0=0$. The IM Tool assumes $E_0=0$ and $P_0=0$ and $V_0=1/\rho_0$, where $\rho_0$ is the initial density in the material database.



## SESAME Tables

SESAME tables are gridded EOS surfaces that are interpolated during problem execution to solve the energy equation. pyKO uses bilinear interpolation on the density-temperature grid.

Supported table formats:
* LANL standard 201+301 format. This table format specficiation is described in {cite:t}`Lyon1992` (<a href="https://github.com/ststewart/aneos-forsterite-2019/blob/master/EOS-docs/Lyon-Johnson-1992-SESAME-database.pdf">PDF file</a>).

```{margin} Table interpolation
Future documentation will discuss pros and cons of different interpolation schemes.
```



## References


```{bibliography}
```

