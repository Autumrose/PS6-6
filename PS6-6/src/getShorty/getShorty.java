package getShorty;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.Set;

import getShorty.getShorty.GraphWeighted.NodeWeighted;

public class getShorty {
	public static void main(String[] args) {
        GraphWeighted graphWeighted = new GraphWeighted(true);
        NodeWeighted zero = new NodeWeighted(0);
        NodeWeighted one = new NodeWeighted(1);
        NodeWeighted two = new NodeWeighted(2);
        NodeWeighted three = new NodeWeighted(3);
        NodeWeighted four = new NodeWeighted(4);
        NodeWeighted five = new NodeWeighted(5);
        NodeWeighted six = new NodeWeighted(6);

        // Our addEdge method automatically adds Nodes as well.
        // The addNode method is only there for unconnected Nodes,
        // if we wish to add any
        graphWeighted.addEdge(zero, one, 8);
        graphWeighted.addEdge(zero, two, 11);
        graphWeighted.addEdge(one, three, 3);
        graphWeighted.addEdge(one, four, 8);
        graphWeighted.addEdge(one, two, 7);
        graphWeighted.addEdge(two, four, 9);
        graphWeighted.addEdge(three, four, 5);
        graphWeighted.addEdge(three, five, 2);
        graphWeighted.addEdge(four, six, 6);
        graphWeighted.addEdge(five, four, 1);
        graphWeighted.addEdge(five, six, 8);

        graphWeighted.DijkstraShortestPath(zero, six);
    }
//	public static void main(String args[]) {
//		Scanner input = new Scanner(System.in);
//		String nextLine = input.nextLine();
//		int intersections = Integer.parseInt(nextLine.substring(0, nextLine.indexOf(' ')));
//		nextLine = nextLine.substring(nextLine.indexOf(' ') + 1);
//		int corridors = Integer.parseInt(nextLine.substring(0, nextLine.length()));
//		GraphWeighted graph = new GraphWeighted(false);
//		int source;
//		int destination;
//		double edge;
//		ArrayList<NodeWeighted> list = new ArrayList<NodeWeighted>();
//		for(int i = 0; i < intersections; i++) {
//			list.add(new NodeWeighted(i));
//		}
//		while (input.hasNext()) {
//			for (int i = 0; i < intersections; i++) {
//				nextLine = input.nextLine();
//				source = Integer.parseInt(nextLine.substring(0, nextLine.indexOf(' ')));
//				nextLine = nextLine.substring(nextLine.indexOf(' ') + 1);
//				destination = Integer.parseInt(nextLine.substring(0, nextLine.length()));
//				nextLine = nextLine.substring(nextLine.indexOf(' ') + 1);
//				edge = Double.parseDouble(nextLine.substring(0, nextLine.length()));
//				graph.addEdge(list.get(source), list.get(destination), edge);
//			}
//			graph.DijkstraShortestPath(list.get(0), list.get(intersections-1));
//			
//			if (input.hasNextLine()) {
//				nextLine = input.nextLine();
//				intersections = Integer.parseInt(nextLine.substring(0, nextLine.indexOf(' ')));
//				nextLine = nextLine.substring(nextLine.indexOf(' ') + 1);
//				corridors = Integer.parseInt(nextLine.substring(0, nextLine.length()));
//			} else
//				return;
//
//		}
//	}

	public static class GraphWeighted {
		private Set<NodeWeighted> nodes;
		private boolean directed;

		GraphWeighted(boolean directed) {
			this.directed = directed;
			nodes = new HashSet<>();
		}

		public void addNode(NodeWeighted... n) {
			// We're using a var arg method so we don't have to call
			// addNode repeatedly
			nodes.addAll(Arrays.asList(n));
		}

		public void addEdge(NodeWeighted source, NodeWeighted destination, double weight) {
			// Since we're using a Set, it will only add the nodes
			// if they don't already exist in our graph
			nodes.add(source);
			nodes.add(destination);

			// We're using addEdgeHelper to make sure we don't have duplicate edges
			addEdgeHelper(source, destination, weight);

			if (!directed && source != destination) {
				addEdgeHelper(destination, source, weight);
			}
		}

		private void addEdgeHelper(NodeWeighted a, NodeWeighted b, double weight) {
			// Go through all the edges and see whether that edge has
			// already been added
			for (EdgeWeighted edge : a.edges) {
				if (edge.source == a && edge.destination == b) {
					// Update the value in case it's a different one now
					edge.weight = weight;
					return;
				}
			}
			// If it hasn't been added already (we haven't returned
			// from the for loop), add the edge
			a.edges.add(new EdgeWeighted(a, b, weight));
		}

		public void printEdges() {
			for (NodeWeighted node : nodes) {
				LinkedList<EdgeWeighted> edges = node.edges;

				if (edges.isEmpty()) {
					System.out.println("Node " + node.n + " has no edges.");
					continue;
				}
				System.out.print("Node " + node.n + " has edges to: ");

				for (EdgeWeighted edge : edges) {
					System.out.print(edge.destination.n + "(" + edge.weight + ") ");
				}
				System.out.println();
			}
		}

		public boolean hasEdge(NodeWeighted source, NodeWeighted destination) {
			LinkedList<EdgeWeighted> edges = source.edges;
			for (EdgeWeighted edge : edges) {
				// Again relying on the fact that all classes share the
				// exact same NodeWeighted object
				if (edge.destination == destination) {
					return true;
				}
			}
			return false;
		}

		public void resetNodesVisited() {
			for (NodeWeighted node : nodes) {
				node.unvisit();
			}
		}

		public void DijkstraShortestPath(NodeWeighted start, NodeWeighted end) {
			// We keep track of which path gives us the shortest path for each node
			// by keeping track how we arrived at a particular node, we effectively
			// keep a "pointer" to the parent node of each node, and we follow that
			// path to the start
			HashMap<NodeWeighted, NodeWeighted> changedAt = new HashMap<>();
			changedAt.put(start, null);

			// Keeps track of the shortest path we've found so far for every node
			HashMap<NodeWeighted, Double> shortestPathMap = new HashMap<>();

			// Setting every node's shortest path weight to positive infinity to start
			// except the starting node, whose shortest path weight is 0
			for (NodeWeighted node : nodes) {
				if (node == start)
					shortestPathMap.put(start, 0.0);
				else
					shortestPathMap.put(node, Double.POSITIVE_INFINITY);
			}

			// Now we go through all the nodes we can go to from the starting node
			// (this keeps the loop a bit simpler)
			for (EdgeWeighted edge : start.edges) {
				shortestPathMap.put(edge.destination, edge.weight);
				changedAt.put(edge.destination, start);
			}

			start.visit();

			// This loop runs as long as there is an unvisited node that we can
			// reach from any of the nodes we could till then
			while (true) {
				NodeWeighted currentNode = closestReachableUnvisited(shortestPathMap);
				// If we haven't reached the end node yet, and there isn't another
				// reachable node the path between start and end doesn't exist
				// (they aren't connected)
				if (currentNode == null) {
					System.out.println("There isn't a path between " + start.n + " and " + end.n);
					return;
				}

				// If the closest non-visited node is our destination, we want to print the path
				if (currentNode == end) {
					System.out.println(
							"The path with the smallest weight between " + start.n + " and " + end.n + " is:");

					NodeWeighted child = end;

					// It makes no sense to use StringBuilder, since
					// repeatedly adding to the beginning of the string
					// defeats the purpose of using StringBuilder
					String path = Integer.toString(end.n);
					while (true) {
						NodeWeighted parent = changedAt.get(child);
						if (parent == null) {
							break;
						}

						// Since our changedAt map keeps track of child -> parent relations
						// in order to print the path we need to add the parent before the child and
						// it's descendants
						path = parent.n + " " + path;
						child = parent;
					}
					System.out.println(path);
					System.out.println("The path costs: " + shortestPathMap.get(end));
					return;
				}
				currentNode.visit();

				// Now we go through all the unvisited nodes our current node has an edge to
				// and check whether its shortest path value is better when going through our
				// current node than whatever we had before
				for (EdgeWeighted edge : currentNode.edges) {
					if (edge.destination.isVisited())
						continue;

					if (shortestPathMap.get(currentNode) + edge.weight < shortestPathMap.get(edge.destination)) {
						shortestPathMap.put(edge.destination, shortestPathMap.get(currentNode) + edge.weight);
						changedAt.put(edge.destination, currentNode);
					}
				}
			}
		}
		public class EdgeWeighted implements Comparable<EdgeWeighted> {

			NodeWeighted source;
			NodeWeighted destination;
			double weight;

			EdgeWeighted(NodeWeighted s, NodeWeighted d, double w) {
				// Note that we are choosing to use the (exact) same objects in the Edge class
				// and in the GraphShow and GraphWeighted classes on purpose - this MIGHT NOT
				// be something you want to do in your own code, but for sake of readability
				// we've decided to go with this option
				source = s;
				destination = d;
				weight = w;
			}

			public String toString() {
				return String.format("(%s -> %s, %f)", source.n, destination.n, weight);
			}

			// We need this method if we want to use PriorityQueues instead of LinkedLists
			// to store our edges, the benefits are discussed later, we'll be using
			// LinkedLists
			// to make things as simple as possible
			public int compareTo(EdgeWeighted otherEdge) {

				// We can't simply use return (int)(this.weight - otherEdge.weight) because
				// this sometimes gives false results
				if (this.weight > otherEdge.weight) {
					return 1;
				} else
					return -1;
			}
		}

		public static class NodeWeighted {
			// The int n and String name are just arbitrary attributes
			// we've chosen for our nodes these attributes can of course
			// be whatever you need
			int n;
			private boolean visited;
			LinkedList<EdgeWeighted> edges;

			NodeWeighted(int n) {
				this.n = n;
				visited = false;
				edges = new LinkedList<>();
			}

			boolean isVisited() {
				return visited;
			}

			void visit() {
				visited = true;
			}

			void unvisit() {
				visited = false;
			}
		}
		private NodeWeighted closestReachableUnvisited(HashMap<NodeWeighted, Double> shortestPathMap) {

			double shortestDistance = Double.POSITIVE_INFINITY;
			NodeWeighted closestReachableNode = null;
			for (NodeWeighted node : nodes) {
				if (node.isVisited())
					continue;

				double currentDistance = shortestPathMap.get(node);
				if (currentDistance == Double.POSITIVE_INFINITY)
					continue;

				if (currentDistance < shortestDistance) {
					shortestDistance = currentDistance;
					closestReachableNode = node;
				}
			}
			return closestReachableNode;
		}
	}
//	static class Graph {
//		int vertices;
//		LinkedList<Edge>[] adjacencylist;
//
//		Graph(int vertices) {
//			this.vertices = vertices;
//			adjacencylist = new LinkedList[vertices];
//			// initialize adjacency lists for all the vertices
//			for (int i = 0; i < vertices; i++) {
//				adjacencylist[i] = new LinkedList<>();
//			}
//		}
//
//		public void addEgde(int source, int destination, double weight) {
//			Edge edge = new Edge(source, destination, weight);
//			adjacencylist[source].addFirst(edge); // for directed graph
//		}
//
//		public void printGraph() {
//			for (int i = 0; i < vertices; i++) {
//				LinkedList<Edge> list = adjacencylist[i];
//				for (int j = 0; j < list.size(); j++) {
//					System.out.println("vertex-" + i + " is connected to " + list.get(j).destination + " with weight "
//							+ list.get(j).weight);
//				}
//			}
//		}
//	}
//
//	static class Edge {
//		int source;
//		int destination;
//		double weight;
//
//		public Edge(int source, int destination, double weight) {
//			this.source = source;
//			this.destination = destination;
//			this.weight = weight;
//		}
//	}
}
