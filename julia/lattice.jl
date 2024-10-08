#%%
using Plots
theme(:juno)

function fibonacci_lattice(n::Int)
    golden_ratio = (1 + sqrt(5)) / 2
    p = golden_ratio - 1 # equals 1 /golden_ratio    
    q = n - 1
    # modf returns (frac, int)
    fibonacci_lattice = [(modf(i / p)[1], i / q) for i in range(0, length=n)]
    return fibonacci_lattice
end
function square_to_disk(point::Tuple{Float64, Float64})
    x, y = point
    theta, r = 2pi * x, sqrt(y)
    return (r * cos(theta), r * sin(theta))
end
function square_to_sphere(point::Tuple{Float64, Float64})
    x, y = point
    theta, phi = 2pi * x, acos(1 - 2y)
    return (cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi))
end
function main()
    n = 1000
    # Define a colormap
    # Create a gradient colormap with n colors
    colors = cgrad(:rainbow, n)
    # Create a grid of points
    lattice = fibonacci_lattice(n)
    # println(lattice)
    plot = scatter(lattice, label="Grid Points", aspect_ratio=:equal, color=colors[1:n])
    # Save the plot to a file (optional)
    savefig(plot, "julia/fibonacci_lattice.png")  # Saves the plot as a PNG file

    lattice_disk = square_to_disk.(lattice)
    plot = scatter(lattice_disk, label="Grid Points", aspect_ratio=:equal, color=colors[1:n])
    savefig(plot, "julia/fibonacci_disk.png")  # Saves the plot as a PNG file

    lattice_sphere = square_to_sphere.(lattice)
    xlims, ylims, zlims = [-1,1], [-1,1], [-1,1]
    plot = scatter3d(
        lattice_sphere, label="Grid Points", aspect_ratio=:equal, color=colors[1:n],
        xlims=xlims, ylims=ylims, zlims=zlims,
        xlabel="X Axis", ylabel="Y Axis", zlabel="Z Axis", 
        legend=false,
        camera=(30, 15)  # Set the azimuthal angle to 30 degrees and elevation to 30 degrees
    )
    savefig(plot, "julia/fibonacci_sphere.png")  # Saves the plot as a PNG file

    return 0
end

# include
# @time main()
if abspath(PROGRAM_FILE) == @__FILE__
    #Test code
    @time main()
end