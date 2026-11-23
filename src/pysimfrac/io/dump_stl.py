import numpy as np
import shutil
import os


def _write_ascii_stl(filename, triangles, name):
    """
    Write triangles to an ASCII STL file.

    Args:
        filename (str): Output file name.
        triangles (np.ndarray): Array of shape (n, 3, 3), the three vertices of
            each triangle. Vertex order sets the facet normal (right-hand rule).
        name (str): Name written on the `solid` line.
    """
    triangles = np.asarray(triangles, dtype=float)
    normals = np.cross(triangles[:, 1] - triangles[:, 0],
                       triangles[:, 2] - triangles[:, 0])
    length = np.linalg.norm(normals, axis=1, keepdims=True)
    # degenerate triangles (e.g. walls where the aperture is zero) get a zero normal
    normals = np.divide(normals, length, out=np.zeros_like(normals),
                        where=length > 0)

    facet = ("  facet normal %.12e %.12e %.12e\n"
             "    outer loop\n"
             "      vertex %.12e %.12e %.12e\n"
             "      vertex %.12e %.12e %.12e\n"
             "      vertex %.12e %.12e %.12e\n"
             "    endloop\n"
             "  endfacet")
    rows = np.hstack([normals, triangles.reshape(-1, 9)])
    with open(filename, 'w') as fp:
        fp.write(f"solid {name}\n")
        np.savetxt(fp, rows, fmt=facet)
        fp.write(f"endsolid {name}\n")


def _grid_triangles(x, y, z, up=True):
    """
    Triangulate a height field, two triangles per grid cell.

    Args:
        x (np.ndarray): Node x coordinates, length nx + 1.
        y (np.ndarray): Node y coordinates, length ny + 1.
        z (np.ndarray): Heights of shape (ny + 1, nx + 1).
        up (bool): If True the normals point in +z, otherwise in -z.

    Returns:
        np.ndarray: Triangles of shape (2 * nx * ny, 3, 3).
    """
    xv, yv = np.meshgrid(x, y)
    p = np.stack([xv, yv, z], axis=-1)
    p00, p10 = p[:-1, :-1], p[:-1, 1:]
    p01, p11 = p[1:, :-1], p[1:, 1:]
    if up:
        tris = [(p00, p10, p11), (p00, p11, p01)]
    else:
        tris = [(p00, p11, p10), (p00, p01, p11)]
    return np.concatenate([np.stack(t, axis=-2).reshape(-1, 3, 3) for t in tris])


def _wall_triangles(xy, z_low, z_high):
    """
    Triangulate a vertical wall between two edges, two triangles per segment.

    Args:
        xy (np.ndarray): (n, 2) horizontal positions along the wall.
        z_low (np.ndarray): Heights of the lower edge at those positions.
        z_high (np.ndarray): Heights of the upper edge at those positions.

    Returns:
        np.ndarray: Triangles of shape (2 * (n - 1), 3, 3). The normals point
        along (direction of travel along `xy`) x (+z).
    """
    low = np.column_stack([xy, z_low])
    high = np.column_stack([xy, z_high])
    l0, l1, h0, h1 = low[:-1], low[1:], high[:-1], high[1:]
    return np.concatenate([np.stack([l0, l1, h1], axis=1),
                           np.stack([l0, h1, h0], axis=1)])


def _side_walls(x, y, z_low, z_high):
    """
    The four vertical walls between two height fields on the same grid,
    with normals pointing out of the domain.

    Args:
        x (np.ndarray): Node x coordinates, length nx + 1.
        y (np.ndarray): Node y coordinates, length ny + 1.
        z_low (np.ndarray or float): Lower heights, (ny + 1, nx + 1) or scalar.
        z_high (np.ndarray or float): Upper heights, (ny + 1, nx + 1) or scalar.

    Returns:
        dict: {'inlet' (y min), 'outlet' (y max), 'left' (x min),
        'right' (x max)} of triangle arrays.
    """
    z_low = np.broadcast_to(z_low, (len(y), len(x)))
    z_high = np.broadcast_to(z_high, (len(y), len(x)))

    def edge(xs, ys):
        xs, ys = np.broadcast_arrays(xs, ys)
        return np.column_stack([xs, ys])

    # each wall is traversed so that (travel direction) x (+z) points outward
    return {
        'inlet': _wall_triangles(edge(x, y[0]), z_low[0, :], z_high[0, :]),
        'outlet': _wall_triangles(edge(x[::-1], y[-1]), z_low[-1, ::-1],
                                  z_high[-1, ::-1]),
        'left': _wall_triangles(edge(x[0], y[::-1]), z_low[::-1, 0],
                                z_high[::-1, 0]),
        'right': _wall_triangles(edge(x[-1], y), z_low[:, -1], z_high[:, -1]),
    }


def dump_stl(self, interior=True, exterior=False, cleanup_run=True):
    """
    Write the fracture surfaces to ASCII STL files.

    Args:
        interior (bool, optional):
            If True, also writes the four side walls of the fracture void
            between the bottom and top surfaces: `inlet.stl` (y = 0),
            `outlet.stl` (y = ly), `left.stl` (x = 0) and `right.stl` (x = lx).
            Defaults to True.

        exterior (bool, optional):
            If True, also writes closed surfaces of the solid blocks on either
            side of the fracture: `bottom_external.stl` (below the bottom
            surface), `top_external.stl` (above the top surface) and
            `whole_block.stl` (both blocks). Defaults to False.

        cleanup_run (bool, optional):
            If True, the files are written into a fresh `stl_files/`
            directory, otherwise into the current directory. Defaults to True.

    Notes:
        `top.stl` and `bottom.stl` are always written. Together with the
        interior walls they form a closed surface around the fracture void,
        with vertices shared along every edge and normals pointing out of the
        void. Each grid cell is split into two triangles.

        The surfaces are periodic (FFT synthesis), with nx samples per period
        in x. The first row/column is repeated at the end so the STL spans
        exactly one period, z(0) == z(lx), on an (nx + 1) x (ny + 1) grid
        with spacing lx / nx = h.

    Example:
        >>> myfrac.dump_stl(interior=True, exterior=True)
    """

    print("--> Dumping fracture surfaces to STL")
    top = np.pad(self.top, ((0, 1), (0, 1)), mode="wrap")
    bottom = np.pad(self.bottom, ((0, 1), (0, 1)), mode="wrap")
    x = np.linspace(0, self.lx, self.nx + 1)
    y = np.linspace(0, self.ly, self.ny + 1)

    out_dir = "."
    if cleanup_run:
        out_dir = "stl_files"
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)

    def write(name, triangles):
        _write_ascii_stl(os.path.join(out_dir, f"{name}.stl"), triangles, name)

    write("top", _grid_triangles(x, y, top, up=True))
    write("bottom", _grid_triangles(x, y, bottom, up=False))

    if interior:
        for name, triangles in _side_walls(x, y, bottom, top).items():
            write(name, triangles)

    if exterior:
        z_min = min(top.min(), bottom.min())
        z_max = max(top.max(), bottom.max())
        # flat outer faces of the blocks, 20% beyond the surfaces
        z_floor = z_min - 0.2 * abs(z_min)
        z_ceiling = z_max + 0.2 * abs(z_max)

        lower = [_grid_triangles(x, y, bottom, up=True),
                 _grid_triangles(x, y, np.full_like(bottom, z_floor), up=False)]
        lower += list(_side_walls(x, y, z_floor, bottom).values())
        lower = np.concatenate(lower)

        upper = [_grid_triangles(x, y, top, up=False),
                 _grid_triangles(x, y, np.full_like(top, z_ceiling), up=True)]
        upper += list(_side_walls(x, y, top, z_ceiling).values())
        upper = np.concatenate(upper)

        write("bottom_external", lower)
        write("top_external", upper)
        write("whole_block", np.concatenate([lower, upper]))

    print("--> Dumping fracture surfaces to STL Complete")
