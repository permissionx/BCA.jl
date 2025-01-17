module Father
export Son
module Son
function f()
    println("THis is son\n")
end
export f
end


end

using .Father
Son.f()