#   figure(
#     diagram(
#       {
#         // let (A, B, C) = ((0,0), (1,0), (2,0))
#         node((0, 0), $i Gamma_a / 2$)
#         edge((-1, -1), (0, 0), "--<|--")
#         edge((-1, 1), (0, 0), "-<|-")
#         edge((0, 0), (1, -1), "-<|-")
#         edge((0, 0), (1, 1), "-<|-")
#       },
#       mark-scale: 1.5,
#     ),
#     caption: [],
#   ),
using Typstry
render(typst"""
    #import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

    #set page(width: auto, height: auto, margin: .5cm)

    #let loop(coord, style, radius: 7mm, ..args) = {
        let θ = 2*calc.atan(2*radius/1pt)
        node(coord, radius: 1pt)
        edge(coord, style, coord, bend: θ, ..args)
    }
    #diagram(
      {
        let (A, B, C) = ((0, 0), (1, 0), (2, 0))
        node(A, $x_1$)
        node(C, $x_2$)
        edge(A, B, "-<|-")
        edge(B, C, "--<|--")
        loop(B, "--|>--")
      },
      mark-scale: 1.5,
    )
    """; output = "cetz_example.pdf");
