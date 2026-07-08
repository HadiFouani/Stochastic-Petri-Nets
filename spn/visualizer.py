import networkx as nx
import matplotlib.pyplot as plt


class PetriNetVisualizer:
    def __init__(self, net):
        self.net = net

    def build_graph(self):
        G = nx.DiGraph()

        for pid, place in enumerate(self.net.places):
            G.add_node(
                f"p{pid}",
                label=place.name,
                kind="place",
            )

        for tid, transition in enumerate(self.net.transitions):
            G.add_node(
                f"t{tid}",
                label=transition.name,
                kind="transition",
            )

            for p, w in transition.pre.items():
                G.add_edge(
                    f"p{p}",
                    f"t{tid}",
                    weight=w,
                )

            for p, w in transition.post.items():
                G.add_edge(
                    f"t{tid}",
                    f"p{p}",
                    weight=w,
                )

        return G

    def draw(self, marking=None, figsize=(12, 8)):
        G = self.build_graph()

        places = [n for n, d in G.nodes(data=True) if d["kind"] == "place"]
        transitions = [n for n, d in G.nodes(data=True) if d["kind"] == "transition"]

        pos = nx.spring_layout(G, seed=42)

        plt.figure(figsize=figsize)

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=places,
            node_shape="o",
            node_size=1800,
            edgecolors="black",
            linewidths=1.5,
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=transitions,
            node_shape="s",
            node_size=1200,
            edgecolors="black",
            linewidths=1.5,
        )
        nx.draw_networkx_edges(
        G,
        pos,
        edgelist=G.edges(),
        arrows=True,
        arrowstyle="-|>",
        arrowsize=22,
        width=1.6,
        connectionstyle="arc3,rad=0.08",
        min_source_margin=20,
        min_target_margin=20,
    )

    
        labels = {}
        for n, d in G.nodes(data=True):
            label = d["label"]

            if d["kind"] == "place" and marking is not None:
                pid = int(n[1:])
                label = f"{label}\n[{marking[pid]}]"

            labels[n] = label

        nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)

        edge_labels = {
            (u, v): d["weight"]
            for u, v, d in G.edges(data=True)
            if d["weight"] != 1
        }

        nx.draw_networkx_edge_labels(
            G,
            pos,
            edge_labels=edge_labels,
            font_size=9,
        )

        plt.axis("off")
        plt.tight_layout()
        plt.show()