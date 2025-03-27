import numpy as np
import cubit
import math
import matplotlib.pyplot as plt

def get_closest_point(rz, r, z):
    """
    Find the closest point in the list to the given point (r, z).
    :param rz: List of points.
    :param r: x-coordinate of the point.
    :param z: y-coordinate of the point.
    :return: Index of the closest point and the distance.
    """
    distances = [
        np.sqrt((rz[i][0] - r) ** 2 + (rz[i][1] - z) ** 2) for i in range(len(rz))
    ]
    return distances.index(min(distances)), min(distances)


def get_grid_points(rz, nr, nz, r_min, r_max, z_min, z_max, max_distance):
    """
    Create a grid of points in the given range.
    :param rz: List of points.
    :param nr: Number of points in the r direction.
    :param nz: Number of points in the z direction.
    :param r_min: Minimum r value.
    :param r_max: Maximum r value.
    :param z_min: Minimum z value.
    :param z_max: Maximum z value.
    :param max_distance: Maximum allowable distance for a grid point.
    :return: List of grid points.
    """
    grid_points = []
    for i in range(nr):
        for j in range(nz):
            r = r_min + i * (r_max - r_min) / (nr - 1)
            z = z_min + j * (z_max - z_min) / (nz - 1)
            ind, dist = get_closest_point(rz, r, z)
            if dist < max_distance:
                grid_points.append([rz[ind][0], rz[ind][1]])
    return np.array(grid_points)


def plot_grid_vs_rz(rz, grid_points):
    """
    Plot the original rz points and the grid points for comparison.
    :param rz: List of [R, Z] coordinates.
    :param grid_points: List of grid points.
    """
    rz = np.array(rz)

    plt.figure(figsize=(10, 8))
    plt.scatter(rz[:, 0], rz[:, 1], s=10, label="Original RZ Points", alpha=0.5)
    plt.scatter(
        grid_points[:, 0],
        grid_points[:, 1],
        color="red",
        s=30,
        label="Grid Points",
        alpha=0.8,
    )
    plt.xlabel("R (m)")
    plt.ylabel("Z (m)")
    plt.title("Comparison of Original RZ Points and Grid Points")
    plt.legend()
    plt.grid()
    plt.axis("equal")
    plt.show()


def mesh_surfaces(
    surfaces_vv, surface_parts, file_name, file_name_output, mesh_folder, grid_output_file, vv_mesh_size
):
    """
    Mesh the surfaces and generate grid points.
    :param surfaces_vv: List of surfaces for the VV group.
    :param surface_parts: Dictionary of surface parts and their mesh sizes.
    :param file_name: Input file name.
    :param file_name_output: Output file name.
    :param mesh_folder: Folder to save mesh files.
    :param grid_output_file: File to save grid points.
    """
    surfaces_vv = list(map(int, surfaces_vv.split()))
    surface_vv = {"surfaces": surfaces_vv, "mesh_size": vv_mesh_size}

    for part in surface_parts:
        surface_parts[part]["surfaces"] = list(
            map(int, surface_parts[part]["surfaces"].split())
        )

    cubit.init([""])
    groups_id = {}

    cubit.cmd(f'open "{file_name}"')
    cubit.cmd("set trimesher geometry sizing off")
    cubit.cmd("set trimesher split overconstrained edges on")

    # Mesh VV surfaces
    surfaces = surface_vv["surfaces"]
    for surface in surfaces:
        cubit.cmd(f"surface {surface} scheme trimesh")
        cubit.cmd(f"surface {surface} size {surface_vv['mesh_size']}")
        cubit.cmd(f"surface {surface} scheme trimesh minimum size {surface_vv['mesh_size']}")
        cubit.cmd(f"mesh surface {surface}")
        cubit.cmd(f'group "vv" add surface {surface}')
    groups_id["vv"] = 2

    # Mesh surface parts
    for pn, part in enumerate(surface_parts):
        surfaces = surface_parts[part]["surfaces"]
        groups_id[part] = pn + 3
        for surface in surfaces:
            cubit.cmd(f"surface {surface} scheme trimesh")
            cubit.cmd(f"surface {surface} size {surface_parts[part]['mesh_size']}")
            cubit.cmd(f"surface {surface} scheme trimesh minimum size {surface_vv['mesh_size']}")
            cubit.cmd(f"mesh surface {surface}")
            cubit.cmd(f'group "{part}" add surface {surface}')

    # Generate mesh and grid points
    nodes = {}
    mesh = {}
    tri_id = 0
    rz = []

    for pn, part in enumerate(surface_parts):
        surfaces = surface_parts[part]["surfaces"]
        mesh[part] = ""
        for surface in surfaces:
            tris = cubit.get_surface_tris(surface)
            for tri in tris:
                tri_id += 1
                n = cubit.get_connectivity("tri", int(tri))
                for i in range(3):
                    if n[i] not in nodes:
                        nodes[n[i]] = cubit.get_nodal_coordinates(int(n[i]))
                mesh[part] += (
                    f"{tri_id} ({nodes[n[0]][0]} {nodes[n[0]][1]} {nodes[n[0]][2]}) "
                    f"({nodes[n[1]][0]} {nodes[n[1]][1]} {nodes[n[1]][2]}) "
                    f"({nodes[n[2]][0]} {nodes[n[2]][1]} {nodes[n[2]][2]}) \n"
                )
        with open(f"{mesh_folder}/{part}.txt", "w") as f:
            f.write(mesh[part])

    # Process VV surfaces
    surfaces = surface_vv["surfaces"]
    mesh["vv"] = ""
    for surface in surfaces:
        tris = cubit.get_surface_tris(surface)
        for tri in tris:
            tri_id += 1
            n = cubit.get_connectivity("tri", int(tri))
            for i in range(3):
                if n[i] not in nodes:
                    nodes[n[i]] = cubit.get_nodal_coordinates(int(n[i]))
                    R = np.sqrt(nodes[n[i]][0] ** 2 + nodes[n[i]][1] ** 2)
                    Z = nodes[n[i]][2]
                    toro = math.degrees(math.atan2(nodes[n[i]][1], nodes[n[i]][0]))
                    if toro < 0:
                        toro += 360
                    rz.append([R, Z, toro * 2 * math.pi / 360])


    rz = np.array(rz)
    cubit.cmd(f'save as "{file_name_output}" overwrite')
    grid_points = get_grid_points(
        rz, 11, 25, min(rz[:, 0]), max(rz[:, 0]), min(rz[:, 1]), max(rz[:, 1]), max_distance=0.1
    )
    plot_grid_vs_rz(rz, grid_points)
    np.savetxt(grid_output_file, grid_points)

if __name__ == "__main__":   
    surfaces_vv = "17254 17252 1684 1685 1686 1687 17253 1657 1658 1659 1660 17255"
    vv_mesh_size = 0.15
    surface_parts = {"limiter-inboard-1": {"surfaces":"17260", "mesh_size": 0.075},
                    "limiter-inboard-2": {"surfaces":"17266", "mesh_size": 0.075},
                    "limiter-inboard-3": {"surfaces":"17272", "mesh_size": 0.075},
                    "limiter-inboard-4": {"surfaces":"17278", "mesh_size": 0.075},
                    "limiter-inboard-5": {"surfaces":"17284", "mesh_size": 0.075},
                    "limiter-inboard-6": {"surfaces":"17290", "mesh_size": 0.075},
                    "limiter-inboard-7": {"surfaces":"17296", "mesh_size": 0.075},
                    "limiter-inboard-8": {"surfaces":"17302", "mesh_size": 0.075},
                    "limiter-inboard-9": {"surfaces":"17308", "mesh_size": 0.075},
                    "limiter-inboard-10": {"surfaces":"17314", "mesh_size": 0.075},
                    "limiter-inboard-11": {"surfaces":"17320", "mesh_size": 0.075},
                    "limiter-inboard-12": {"surfaces":"17326", "mesh_size": 0.075},
                    "limiter-inboard-13": {"surfaces":"17332", "mesh_size": 0.075},
                    "limiter-inboard-14": {"surfaces":"17338", "mesh_size": 0.075},
                    "limiter-inboard-15": {"surfaces":"17344", "mesh_size": 0.075},
                    "limiter-inboard-16": {"surfaces":"17350", "mesh_size": 0.075},
                    "limiter-inboard-17": {"surfaces":"17356", "mesh_size": 0.075},
                    "limiter-inboard-18": {"surfaces":"17362", "mesh_size": 0.075},
                    "limiter-outboard-1": {"surfaces":"17368", "mesh_size": 0.075},
                    "limiter-outboard-2": {"surfaces":"17374", "mesh_size": 0.075},
                    "limiter-outboard-3": {"surfaces":"17380", "mesh_size": 0.075},
                    "limiter-outboard-4": {"surfaces":"17386", "mesh_size": 0.075},
                    "limiter-outboard-5": {"surfaces":"17392", "mesh_size": 0.075},
                    "limiter-outboard-6": {"surfaces":"17398", "mesh_size": 0.075},
                    "limiter-outboard-7": {"surfaces":"17404", "mesh_size": 0.075},
                    "limiter-outboard-8": {"surfaces":"17410", "mesh_size": 0.075},
                    "limiter-outboard-9": {"surfaces":"17416", "mesh_size": 0.075},
                    "limiter-outboard-10": {"surfaces":"17422", "mesh_size": 0.075},
                    "limiter-outboard-11": {"surfaces":"17428", "mesh_size": 0.075},
                    "limiter-outboard-12": {"surfaces":"17434", "mesh_size": 0.075},
                    "limiter-outboard-13": {"surfaces":"17440", "mesh_size": 0.075},
                    "limiter-outboard-14": {"surfaces":"17446", "mesh_size": 0.075},
                    "limiter-outboard-15": {"surfaces":"17452", "mesh_size": 0.075},
                    "limiter-outboard-16": {"surfaces":"17458", "mesh_size": 0.075},
                    "limiter-outboard-17": {"surfaces":"17464", "mesh_size": 0.075},
                    "limiter-outboard-18": {"surfaces":"17470", "mesh_size": 0.075}}

    file_name = "first_wall_tracer.cub" # input geometry 
    file_name_output = "first_wall_tracer_out.cub" # output geometry
    mesh_folder = "mesh" # folder to save mesh files
    grid_output_file = 'grid_points.txt' # file to save grid points

    mesh_surfaces(surfaces_vv, surface_parts, file_name, file_name_output, mesh_folder, grid_output_file, vv_mesh_size)
