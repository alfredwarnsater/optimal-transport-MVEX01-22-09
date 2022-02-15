using Plots

a = [1, 2, 3, 2, 1]
b = [1, 2, 3, 4, 5]

print(a)
p1 = bar(a, orientation = :horizontal, xflip = true)
p2 = bar(b, orientation = :horizontal, xflip = false)


plot(p1, p2, layout=(2,1))