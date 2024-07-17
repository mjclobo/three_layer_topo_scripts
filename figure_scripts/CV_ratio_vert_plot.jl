# loading modules
using PyPlot, CSV, DataFrames, Printf

using JLD2, FileIO

# controls ratio of interfaces' stratification
gammas = round.(collect(range(0.2,3,15)),sigdigits=2) #

# controls ratio of interface shears
alphas = round.(collect(range(1,5,21)),sigdigits=2)

# topo parameters
# h0s_ext = [-10 .^reverse(collect(range(-7,-3,6))); 0. ; 10 .^collect(range(-7,-3,6))]
# h0s_fill = [ - collect(range(10^-4.5,10^-3.2,10)); collect(range(10^-4.5,10^-3.2,10))]

# h0s_ext = sort([h0s_ext; h0s_fill]);

h0s_ext = collect(range(-10^-3,10^-3,30));
# h0s_ext = [h0s_ext[1:15]; 0.; h0s_ext[16:end]];


kts = [0]  # zero wavenumber for sloping topography

Ng = length(gammas); Na = length(alphas); Nh_ext = length(h0s_ext); Nk = length(kts);

H = [500.,1000.,2500.]; Nz = length(H);
H_grid = -[0, H[1], sum(H[1:2]), sum(H)];

f0 = 8.3e-5

# get file names
dat_path = "../data/data_jld2_02/"
# dat_path = "../data/data_512_01/"

all_files = readdir(dat_path)

Lx = 1000.e3

filename(gam,alp,h0,k) = "threelayer_power_iter_gamma"*string(gam)*"_alpha"*string(alp)*"_h0"*string(round(h0*Lx,digits=9))*"_kt"*string(Int(k))*"_res256.jld2"

###########################################################################

# plot VF32 as function of h0 and alpha for a given gamma; Charney-favoring
VF32_gam3 = zeros(Nh_ext,Na)
VF52_gam3 = zeros(Nh_ext,Na)
ED_gam3 = zeros(Nh_ext,Na)
FCIC_gam3 = zeros(Nh_ext,Na)

gam = gammas[3];
kt = kts[1];

for (i,h0) in enumerate(h0s_ext)
    for (j,alp) in enumerate(alphas)
        fil = findall(x->x==filename(gam,alp,h0,kt),all_files)
        if isempty(fil)==false
            dat_in = load(dat_path*filename(gam,alp,h0,kt));
            VF32_gam3[i,j] = dat_in["csv_data"]["VF32"][end]
            VF52_gam3[i,j] = dat_in["csv_data"]["VF52"][end]
            ED_gam3[i,j] = dat_in["csv_data"]["Ekman_drag"][end]

            h2oh3 = dat_in["csv_data"]["H"][2] / dat_in["csv_data"]["H"][3]
            S52 = f0* dat_in["csv_data"]["U"][2]/((9.81/1000.0)*(dat_in["csv_data"]["rho"][3] - dat_in["csv_data"]["rho"][2]))

            if (alp/gam) > (1 + h2oh3 * (1 - h0/S52))
                FCIC_gam3[i,j] = 1.0
            else
                 FCIC_gam3[i,j] = NaN
            end
        else
            println("h: "*string(h0))
            println("alpha: "*string(alp))
        end
    end
end

# plot VF32 as function of h0 and alpha for a given gamma; Eady-favoring
VF32_gam13 = zeros(Nh_ext,Na)
VF52_gam13 = zeros(Nh_ext,Na)
ED_gam13  = zeros(Nh_ext,Na)
FCIC_gam13  = zeros(Nh_ext,Na)

gam = gammas[13];
kt = kts[1];

for (i,h0) in enumerate(h0s_ext)
    for (j,alp) in enumerate(alphas)
        fil = findall(x->x==filename(gam,alp,h0,kt),all_files)
        if isempty(fil)==false
            dat_in = load(dat_path*filename(gam,alp,h0,kt));
            VF32_gam13[i,j] = dat_in["csv_data"]["VF32"][end]
            VF52_gam13[i,j] = dat_in["csv_data"]["VF52"][end]
            ED_gam13[i,j] = dat_in["csv_data"]["Ekman_drag"][end]

            h2oh3 = dat_in["csv_data"]["H"][2] / dat_in["csv_data"]["H"][3]
            S52 = f0* dat_in["csv_data"]["U"][2]/((9.81/1000.0)*(dat_in["csv_data"]["rho"][3] - dat_in["csv_data"]["rho"][2]))

            if (alp/gam) >= 1
                if (alp/gam) > (1 + h2oh3 * (1 - h0/S52))
                    FCIC_gam13[i,j] = 1.0
                else
                     FCIC_gam13[i,j] = NaN
                end
            else
                 FCIC_gam13[i,j] = NaN
            end
        else
            println("h: "*string(h0))
            println("alpha: "*string(alp))
        end
    end
end

x_mid = Int(floor(length(h0s_ext)/2)+1)

hxtix=[L"-10e".*string.(log10.(abs.(h0s_ext[1:x_mid-1]))) ; L"0.0" ; L"10e".*string.(log10.(abs.(h0s_ext[x_mid+1:end]))) ]

nlayers = 3                 # number of layers
f0 = f₀= 8.3e-5  # 1.0e-4            # Coriolis param [s^-1]
g = 9.81                    # gravity
H = [500., 1000., 2500.]     # the rest depths of each layer; [250., 500., 3000.]

H32 = 0.5 * (H[1] + H[2])
H52 = 0.5 * (H[2] + H[3])

U = [0.05,0.0255,0.0]

# setting base density profile
rho = ρ = [0.0, 1025.0, 1025.75]         # the density of each layer

U1 = U[2] .+ alphas*(H32/H52)*(U[2]-U[3])
# rho1 = rho[2] - (abs(U[2]-U[1])*(rho[3]-rho[2]))/abs(U[2]-U[3])/gamma
rho1g3 = rho[2] - gammas[3] * (rho[3] - rho[2]) * (H32/H52)

rho1g13 = rho[2] - gammas[13] * (rho[3] - rho[2]) * (H32/H52)

S52g3 = -(f0*rho1g3/g)*((U[2] .- U[3])./(rho[2]-rho[3]))

S52g13 =  -(f0*rho1g13/g)*((U[2] .- U[3])./(rho[2]-rho[3]))

# log10.(abs.(S52g13./h0s_ext))

###########################################################################
matplotlib[:rc]("text", usetex= "true")
matplotlib[:rc]("text.latex", preamble= "\\usepackage{amsmath,amsthm}")

PyPlot.matplotlib[:rc]("hatch",linewidth=1.5)
PyPlot.matplotlib[:rc]("hatch",color="red")


using PlotUtils

fig,ax = PyPlot.subplots(2,1,figsize=(3.0,5.0))
fig.tight_layout(pad=0.0)
ax1=ax[1]; ax2=ax[2];

f_size = 18.
l_size = 12.

u_lim = maximum(log10.((VF32_gam3 ./ VF52_gam3)))
l_lim = minimum(log10.((VF32_gam13 ./ VF52_gam13))) # -u_lim

# ./ (alphas'/gammas[3])
# pc1 = ax1.pcolormesh(h0s_ext/S52g3,alphas/gammas[3],log10.((VF32_gam3 ./ VF52_gam3))', cmap="Spectral", vmin=l_lim, vmax=u_lim)
pc1 = ax1.contour(h0s_ext/S52g3,alphas/gammas[3],log10.((VF32_gam3 ./ VF52_gam3))', vmin=l_lim, vmax=u_lim, colors="black",linewidths=1.5)
ax1.pcolor(h0s_ext/S52g13,alphas/gammas[3], FCIC_gam3',hatch="/", alpha=0.)
ax1.tick_params(labelsize=l_size)
ax1.set_ylim([alphas[1]/gammas[3],alphas[end]/gammas[3]])
ax1.set_xlim(h0s_ext[1]/S52g3, h0s_ext[end]/S52g3)
ax1.text(0.07, 0.925, L"\mathbf{(a)} ", fontsize=19, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes);  # pretty off-white "#F5F5DC"
ax1.text(0.2, 0.8, L"\gamma=0.6", fontsize=12, color="red", horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes, bbox=Dict("facecolor" => "white", "edgecolor" => "none"));
cl1 = ax1.clabel(pc1, pc1.levels, inline=true, fontsize=10, manual=[(-2.5,2),(-1.75,2.1),(-1,3),(-1,4.25),(2,4),(2,6),(2.1,7.5)])
ax1.set_xticklabels([])
[txt.set_backgroundcolor("white") for txt in cl1]
# [txt.set_bbox(Dict(facecolor => "white", "edgecolor" => "none", pad=0)) for txt in cl1]

# ./ (alphas'/gammas[13])
#pc2 = ax2.pcolormesh(h0s_ext/S52g13,alphas/gammas[13],log10.((VF32_gam13 ./ VF52_gam13))',cmap="Spectral", vmin=l_lim, vmax=u_lim)
pc2 = ax2.contour(h0s_ext/S52g13,alphas/gammas[13],log10.((VF32_gam13 ./ VF52_gam13) ./ (alphas'/gammas[13]))', vmin=l_lim, vmax=u_lim,colors="black",linewidths=1.5)
ax2.pcolor(h0s_ext/S52g13,alphas/gammas[13], FCIC_gam13', hatch="/", alpha=0.)
ax2.tick_params(labelsize=l_size)
ax2.set_ylim([alphas[1]/gammas[13],alphas[end]/gammas[13]])
ax2.set_xlim(h0s_ext[1]/S52g13, h0s_ext[end]/S52g13)
ax2.text(0.08, 0.925, L"\mathbf{(b)}", fontsize=19, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);
ax2.text(0.2, 0.8, L"\gamma = 2.6", color="red", fontsize=12, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);
ax2.plot([h0s_ext[1]/S52g13, h0s_ext[end]/S52g13],[1,1],color="cyan",linestyle="dashed",linewidth=2.)
cl2 = ax2.clabel(pc2, pc2.levels, inline=true, fontsize=10, manual=[(-2.9,0.5),(-2,0.7),(-1,0.75),(0.,0.98),(-2,1.2),(0,1.5),(1.75,1.5),(2.8,1.5)])
[txt.set_backgroundcolor("white") for txt in cl2]
# [txt.set_bbox(Dict(facecolor => "white", "edgecolor" => "none", pad=0)) for txt in cl2]

#=
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(pc2, cax=cbar_ax)
cb.ax.tick_params(labelsize=l_size)
cb.ax.yaxis.offsetText.set(size=l_size)
cb.ax.set_yticklabels([latexstring(string.(round.(i; digits=2))) for i in cb.get_ticks()])
=#


fig.text(-0.1, 0.5, L"S_{3/2} / S_{5/2}", va="center", rotation="vertical", fontsize=f_size)
fig.text(0.5, -0.03, L"\partial_y \eta_b / S_{5/2}", va="center", ha="center", rotation="horizontal", fontsize=f_size)
fig.text(0.5, 1.05, L"\mathrm{log}_{10} \left ( \frac{C_{\mathrm{V,}3/2}}{C_{\mathrm{V,}5/2}}  \right )", va="center", ha="center", rotation="horizontal", fontsize=f_size)

savefig("./CV_ratio_vert.png",bbox_inches="tight")

savefig("./CV_ratio_vert.pdf",bbox_inches="tight")

































