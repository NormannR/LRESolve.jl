# LRESolve
Solving Systems of Linear Rational Expectations Equations in Julia

## Installation

## Methods

### Sims (2001)
[Sims (2001)](https://ideas.repec.org/c/dge/qmrbcd/11.html) solves LRE systems of the form

![image](http://www.sciweavers.org/tex2img.php?eq=%5CGamma_0%20x_%7Bt%7D%20%3D%20C%20%2B%20%5CGamma_1%20x_%7Bt-1%7D%20%2B%20%5CPsi%20z_t%20%2B%20%5CPi%20%5Ceta_t&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0)

where 

- ![alt text](http://www.sciweavers.org/tex2img.php?eq=x&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0) is the vector of endogenous variables
- ![alt text](http://www.sciweavers.org/tex2img.php?eq=z&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0) is the vector of exogenous shocks
- η is the vector of expectation errors

The solution verifies

![image](http://www.sciweavers.org/tex2img.php?eq=x_t%20%3D%20%5CTheta%20%2B%20%5CTheta_0%20z_t%20%2B%20%5CTheta_1%20x_%7Bt-1%7D&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0[/img])

To solve a LRE system using this method
1. Define the model through the `ModelSims` structure. The syntax is typically

```julia
M = ModelSims(Γ₀,Γ₁,C,Ψ,Π)
```

2. Call the `solve_sims` method over the newly created model
```julia
Θ, Θ₀, Θ₁ = solve_sims(M)
```
