function test1()
n = 1000000
a = rand(Float32, n)
b = rand(Float64, n)
c = rand(Float16, n)
a2 = sort(a)
b2 = sort(b)
c2 = sort(c)

res1 = 0.0
res2 = 0.0
res3 = 0.0
res4 = 0.0
res5 = 0.0
res6 = 0.0

for i in 1:n
    res1 = res1 + a[i]
    res2 = res2 + b[i]
    res3 = res3 + a2[i]
    res4 = res4 + b2[i]
    res5 = res5 + c[i]
    res6 = res6 + c2[i]
end
println("a  : $res1, type : $(typeof(a))")
println("a2 : $res3, type : $(typeof(a2))")
println("b  : $res2, type : $(typeof(b))")
println("b2 : $res4, type : $(typeof(b2))")
println("c  : $res5, type : $(typeof(c))")
println("c2 : $res6, type : $(typeof(c2))")
end

function test2()
n = 1000000
a = rand(Float128, n)
b = convert(Array{Float32, 1}, a)
c = convert(Array{Float64, 1}, a)
a2 = sort(a)
b2 = sort(b)
c2 = sort(c)

res1::Float128=0.0
res2::Float128=0.0
res3::Float32=0.0
res4::Float32=0.0
res5::Float64=0.0
res6::Float64=0.0

for i in 1:n
    res1 = res1 + a[i]
    res2 = res2 + a2[i]
    res3 = res3 + b[i]
    res4 = res4 + b2[i]
    res5 = res5 + c[i]
    res6 = res6 + c2[i]
end
println("a  : $res1, type : $(typeof(res1))")
println("a2 : $res2, type : $(typeof(res2))")
println("b  : $res3, type : $(typeof(res3))")
println("b2 : $res4, type : $(typeof(res4))")
println("c  : $res5, type : $(typeof(res5))")
println("c2 : $res6, type : $(typeof(res6))")
end

function test3()
n = 1000000
a = rand(0.0:1e-16:1.0, n)
a2 = sort(a)

res1 = 0.0
res2 = 0.0

for i in 1:n
    res1 = res1 + a[i]
    res2 = res2 + a2[i]
end
println("a  : $res1, type : $(typeof(res1))")
println("a2 : $res2, type : $(typeof(res2))")
end
