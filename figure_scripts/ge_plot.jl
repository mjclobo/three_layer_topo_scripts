# growth rate in gamma-eta space

# load modules
include("../LinStab/mjcl_stab.jl")
using .LinStab

using PyPlot, Measures, LaTeXStrings, JLD2, FileIO, Peaks

matplotlib[:rcParams]["axes.unicode_minus"]=false

PyPlot.matplotlib[:rc]("text", usetex=true) 
PyPlot.matplotlib[:rc]("text.latex", preamble= "\\usepackage{amsmath,amsthm}")


# set static vars
L = Lx = Ly = 1000.e3                   # domain size [m]
beta = β = 0 # 1.14052e-11              # the y-gradient of planetary PV
Ny = Nx = n = 256               # 2D resolution = n²

nlayers = 3                 # number of layers
f0 = f₀= 8.3e-5            # Coriolis param [s^-1] 
g = 9.81                    # gravity
H = [500., 1000., 2500.]     # the rest depths of each layer; [250., 500., 3000.] 
Nz = length(H)

U = [0.05,0.025,0.0]

# controls ratio of interface shears
alpha = 2.8


# set gamma range
gammas = [4.0,1.9,0.89];


Ng = length(gammas)

# setting density profile as function of gamma
rho = ρ = [1024.0, 1025.0, 1025.75]         # the density of each layer

V = zeros(nlayers)

# Linear bottom drag
μ = 0 # f0*10/2/H[3]            # d_E is a typical value of 5-10 m in the ocean [Weatherly and Martin, 1978; Perlin et al., 2007]

# step thru gamma while storing vertical structure of most unstable mode
Nk = round(Int,Nx/2)
sigma_gamma = zeros((Ng,Nk))

h0s = collect(range(-10^-3,10^-3,101))

Nh = length(h0s)

Qy_gamma_alpha  = zeros((Nh,Ng,Nz))

H32 = 0.5 * (H[1] + H[2])
H52 = 0.5 * (H[2] + H[3])


S32 =zeros(Ng)



U[1] = U[2] + alpha*(H32/H52)*(U[2]-U[3])

# rho1 = rho[2] - (abs(U[2]-U[1])*(rho[3]-rho[2]))/abs(U[2]-U[3])/gamma

rho1g3 = rho[2] - gammas[1] * (rho[3] - rho[2]) * (H32/H52)

rho1g13 = rho[2] - gammas[1] * (rho[3] - rho[2]) * (H32/H52)


S52g3 = zeros((Nh,Ng))


S52g13 =  -(f0*rho1g13/g)*((U[2] .- U[3])./(rho[2]-rho[3]))



###########################################################
a = load("./LSA_data/ge_sigma_gamma_eta.jld")
max_k = a["sigma_gamma_eta"]

##########################################################
down_samp = 5

maxes = zeros(size(max_k))
for i in range(1,size(max_k)[1])
    for j in range(1,size(max_k)[2])
        if mod(i,down_samp) == 0
            s = max_k[i,j,:]
            q = argmaxima(s)
            @. maxes[i,j,q] = h0s[i]/S52g13
            @. maxes[i,j,maxes[i,j,:] == 0] = NaN
        else
            @. maxes[i,j,:] = NaN
        end
    end
end

###########################################################

using SkipNan

function log_zeros(vec)
    vec[vec .<= (0.01 / 3600 / 24 / 7)] .= NaN
    a = vec .* 3600 .* 24 .* 7;
    # a = log.(vec)
    # a[a .== -Inf] .= 0.0
    return a
end

# growth rates as function of h_y
eta = 0
eve1,eva1,max_eve1,max_eve_phase1,max_eva1,k_x,k_y,qx1,qy1,rd1 = LinStab.lin_stab(U,V,H,beta,eta,Nx,Ny,rho,f0,g,Float64(Lx),Float64(Ly))
k_xr = k_x[round(Int,Nx/2+1):end]

rd40 = 36307
rd19 = 33115
rd89 = 32010

###########################################################

f_size = 36.
l_size = 28.

fig,ax = PyPlot.subplots(1,3,figsize=(16,5),width_ratios=[1.0,1.0,1.3])
fig.tight_layout(pad=3.0)
ax1=ax[1]; ax2=ax[2]; ax3=ax[3];

x_max = 2.5;

ind1 = 1;
ind2 = 2;
ind3 = 3;

max_sig = maximum(skipnan(log_zeros(max_k[:,ind3,1:32])))
ind = max_k[:,2,1:32] .> 0.0;
min_sig = minimum(skipnan(log_zeros(max_k[:,ind1,1:32])))


pc1=ax1.pcolormesh(k_xr*rd40,h0s/S52g13,log_zeros(max_k[:,ind1,:]),vmin=min_sig,vmax=max_sig,cmap=PyPlot.cm.cividis,linewidth=0,rasterized=true)
for i in range(1,size(max_k)[1])
    ax1.plot(k_xr*rd40,maxes[i,ind1,:],"ro",markersize=4.0)
end
ax1.tick_params(labelsize=10.)
ax1.set_ylabel(L"\partial_y \eta_b / S_{5/2}", fontsize=f_size)
ax1.set_title(L"\sigma { \scriptstyle [\mathrm{week}^{-1}] } \ \ {\footnotesize \left (\frac{S_{3/2} }{ S_{5/2}} = 0.7 \right ) }", fontsize = f_size, y=1.1)
ax1.tick_params(axis="both", which="major", labelsize=l_size)
ax1.set_xlim(0,x_max)
ax1.text(0.1, 0.9, L"\mathbf{(a)}", fontsize=42, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes);


pc2=ax2.pcolormesh(k_xr*rd19,h0s/S52g13,log_zeros(max_k[:,ind2,:]),vmin=min_sig,vmax=max_sig,cmap=PyPlot.cm.cividis,linewidth=0,rasterized=true)
for i in range(1,size(max_k)[1])
    ax2.plot(k_xr*rd19,maxes[i,ind2,:],"ro",markersize=4.0)
end
ax2.tick_params(labelsize=10.)
ax2.set_yticks([])
ax2.set_xlabel(L"k_x L_{\mathrm{d,}1} ", fontsize=f_size)
ax2.set_title(L"\sigma { \scriptstyle [\mathrm{week}^{-1}] } \ \ {\footnotesize \left (\frac{S_{3/2} }{ S_{5/2}}  = 1.5 \right )}", fontsize = f_size, y=1.1)
ax2.tick_params(axis="both", which="major", labelsize=l_size)
ax2.set_xlim(0, x_max)
ax2.text(0.1, 0.9, L"\mathbf{(b)}", fontsize=42, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);


pc3=ax3.pcolormesh(k_xr*rd89,h0s/S52g13,log_zeros(max_k[:,ind3,:]),vmin=min_sig,vmax=max_sig,cmap=PyPlot.cm.cividis,linewidth=0,rasterized=true)
for i in range(1,size(max_k)[1])
    ax3.plot(k_xr*rd89,maxes[i,ind3,:],"ro",markersize=4.0)
end
ax3.tick_params(labelsize=10.)
ax3.set_yticks([])
ax3.set_title(L"\sigma { \scriptstyle [\mathrm{week}^{-1}] } \ \ {\footnotesize \left (\frac{S_{3/2} }{ S_{5/2}}  = 3.1 \right )}", fontsize = f_size, y=1.1)
ax3.tick_params(axis="both", which="major", labelsize=l_size)
ax3.set_xlim(0, x_max)
cb = fig.colorbar(pc3, shrink=0.85, pad=0.1)
cb.ax.tick_params(labelsize=l_size)
# cb.ax.yaxis.offsetText.set(size=l_size)
# cb.ax.yaxis.set_offset_position("left")                         
ax3.text(0.1, 0.9, L"\mathbf{(c)}", fontsize=42, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);

pc1.set_edgecolor("face")
pc2.set_edgecolor("face")
pc3.set_edgecolor("face")

savefig("./figs/ge_wave_growth.png",bbox_inches="tight")


savefig("./figs/ge_wave_growth.pdf",bbox_inches="tight")


###########################################################

