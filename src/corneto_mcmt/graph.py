import corneto


def example_graph() -> corneto.Graph:
    g = corneto.Graph()
    g.add_edge("A", "B")
    g.add_edge("B", "C")
    return g
