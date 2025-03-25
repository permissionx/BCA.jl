function Dump(simulator::Simulator, filename::String, step::Int64, isAppend::Bool=false)
    if !simulator.isOrthogonal        
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    write_flag = isAppend ? "a" : "w"
    open(filename, write_flag) do file
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(simulator.numberOfAtoms), "\n")
        write(file, "ITEM: BOX BOUNDS ")
        for d in 1:3
            if simulator.periodic[d] 
                write(file, "pp ")
            else
                write(file, "ff ")
            end
        end
        write(file, "\n")
        for d in 1:3
            write(file, "0 $(simulator.box.vectors[d,d])\n")
        end
        write(file, "ITEM: ATOMS id type x y z vx vy vz energy cx cy cz\n")
        for atom in simulator.atoms
            if atom.isAlive
                write(file, "$(atom.index) $(atom.type) \
                $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
                $(atom.velocityDirection[1]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[2]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[3]*sqrt(2*atom.mass*atom.energy)) \
                $(atom.energy) \
                $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3])\n")
            end
        end 
    end
end

