# Model description

I am modeling reflection (red) using the following, where $x$ is reflection and 
$y$ is the satellite observation. 

$$
logit(x_r) = \alpha + \beta\psi_r
$$

$$
y_r \sim N(x_r, \sigma^2)
$$


Here, $r'$ is the set of all upstream observations, and $w_s$ is the weight for 
upstream site $s$, derived from flow rate.

$\psi_r$ represents spatial variation and is modeled using the Ornstein-Uhlenbeck 
process, where correlation and variance is derived from a weighted average of 
upstream sites.

$$
\psi_r = \sum_{s\in r'}w_{sr}\rho_{sr}\psi_s + \epsilon_r
$$

$$
\epsilon_r \sim Normal(0, \sigma^2_r) 
$$

$\sigma^2_r$ is the spatial variance for each reach $r$, given by

$$
\sigma^2_r = \sum_{s\in r`}w_s (1-e^{-2\theta_v d_s})
$$

where $d_s$ is the distance from observation $r$ to upstream observation $s$.

$\rho_s$ is the correlation due to spatial similarity:

$$
\rho_s = e^{-\theta_v d_s}.
$$

We can express the joint distribution of $\psi$ as a multivariate normal distribution 
in terms of a path matrix $\Gamma$ and individual variance component, with diagonal 
elements corresponding to $\sigma^2_r$ defined above: 

$$
\boldsymbol \psi \sim MVN(0, \Sigma) \\
\Sigma  = (I - \Gamma)^{-1}V(I-\Gamma^\top)
$$

## With covariates

We can also add covariates to the spatial process as 

$$
\psi_r = \sum_{s\in r'}w_{sr}\rho_{sr}\psi_s + \beta_c C_r + \epsilon_r
$$

where $C_r$ is the vector of covariates for $\psi_r$ and $\beta$ is the vector of 
coefficients. Then the vectorized form is 

$$
\boldsymbol \psi = \Gamma \boldsymbol \psi + \boldsymbol \beta_c \boldsymbol C +\boldsymbol \epsilon 
$$

and solving for $\boldsymbol \psi$, we get

$$
\boldsymbol \psi - (\boldsymbol I - \boldsymbol \Gamma)^{-1}\boldsymbol \beta \boldsymbol C  \sim (\boldsymbol I - \boldsymbol \Gamma)^{-1}MVN(0, \boldsymbol \Sigma)
$$

$$
\boldsymbol \Sigma  = (\boldsymbol I - \boldsymbol \Gamma)^{-1}\boldsymbol V((\boldsymbol I- \boldsymbol \Gamma)^\top)^{-1}
$$



## Add temporal component

Adding conditional spatio-temporal autocorrelation.
We add a spatio-temporal random variable, $\omega$. 
Let $\boldsymbol{\psi} \sim MVN(0, \Sigma)$. Then, 

$$
logit(\boldsymbol{x}_t) = \alpha + \beta_\psi\boldsymbol{\psi} + \beta_\omega\boldsymbol{\omega}_t
$$

where $\boldsymbol \psi$ is as defined above. The spatio-temporal component, 
$\omega$, has the same spatial correlation as $\psi$ and an AR-1 temporal 
autocorrelation structure: 

$$
\omega_t \sim     \begin{cases}
        MVN(0, \Sigma) & \text{if } t = 1\\
        MVN(\rho_\omega \omega_{t-1}, \Sigma) & \text{if } t > 1
    \end{cases}
$$

Observations are modeled

$$
\boldsymbol y_t \sim N(\boldsymbol x_r, sigma^2)
$$

