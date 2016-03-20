import numpy as np
import struct
import sys
import yt

def write_cube_ndfield(data, filename):
    print "Writing grid with dims:", data.ActiveDimensions
    with file(filename, "wb") as f:
        # dummy int(4B)
        f.write(struct.pack(">i", 0))
        # tag char(1B) x 16
        f.write(struct.pack("16s", "NDFIELD"))
        # dummy int(4B)
        f.write(struct.pack(">i", 0))
        # dummy int(4B)
        f.write(struct.pack(">i", 0))
        # some fucking comment char(1B) x 80
        f.write(struct.pack("80s", "seriously"))
        # ndims int(4B)
        f.write(struct.pack(">i", 3))
        # dims int(4B) x 20
        for i in range(3):
            f.write(struct.pack(">i", data.ActiveDimensions[i]))
        for i in range(17):
            f.write(struct.pack(">i", 1))
        # fdims_index int(4B)
        f.write(struct.pack(">i", 0))
        # datatype int(4B)
        f.write(struct.pack(">i", 512)) # double
        # x0 double(8B) x 20
        for i in range(3):
            f.write(struct.pack(">d", data.ds.domain_left_edge[i].in_units("code_length")))
        for i in range(17):
            f.write(struct.pack(">d", 0.))
        # delta double(8B) x 20
        for i in range(3):
            f.write(struct.pack(">d", data.ds.domain_right_edge[i].in_units("code_length")))
        for i in range(17):
            f.write(struct.pack(">d", 0.))
        # dummy ext char(1B) x 160
        f.write(struct.pack("160s", "whatever"))
        # dummy int(4B)
        f.write(struct.pack(">i", 0))
        # dummy int(4B)
        f.write(struct.pack(">i", 0))
        # data size of datatype x N (double for us)
        field_data = cg["density"] + cg["dark_matter_density"]
        pbar = yt.get_pbar("Writing ndfield file:", np.prod(data.ActiveDimensions))
        q = 0
        for i in range(data.ActiveDimensions[0]):
            for j in range(data.ActiveDimensions[1]):
                for k in range(data.ActiveDimensions[2]):
                    f.write(struct.pack(">d", field_data[k][j][i]))
                    q += 1
                    pbar.update(q)
        pbar.finish()
        # dummy int(4B)
        f.write(struct.pack(">i", 0))
        f.close()

if __name__ == "__main__":
    ds_fn = sys.argv[1]
    output_fn = sys.argv[2]
     
    ds = yt.load(ds_fn)
    cg_level = 0
    cg = ds.covering_grid(cg_level, ds.domain_left_edge, 2**cg_level * ds.domain_dimensions)

    write_cube_ndfield(cg, output_fn)
