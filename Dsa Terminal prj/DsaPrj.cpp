
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <queue>
#include <stack>
#include <limits>

using namespace std;

#define MAX_VERTICES 1000
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct City
{
    string cityName;
    double latitude;
    double longitude;
    string countryName;
    int population;
    int neighbors[MAX_VERTICES];

    bool operator==(const City &other) const
    {
        return cityName == other.cityName && countryName == other.countryName &&
               latitude == other.latitude && longitude == other.longitude &&
               population == other.population;
    }
};
class HashTable
{
private:
    static const int TABLE_SIZE = 1000;
    struct Node
    {
        string key;
        int value;
        Node *next;
    };

    Node *table[TABLE_SIZE];

    int hash(const string &key)
    {
        int hashValue = 0;
        for (char c : key)
        {
            hashValue = (hashValue * 31 + static_cast<int>(c)) % TABLE_SIZE;
        }
        return hashValue;
    }

public:
    HashTable()
    {
        for (int i = 0; i < TABLE_SIZE; i++)
        {
            table[i] = nullptr;
        }
    }

    void insert(const string &key, int value)
    {
        int index = hash(key);

        Node *newNode = new Node;
        newNode->key = key;
        newNode->value = value;
        newNode->next = nullptr;

        if (table[index] == nullptr)
        {
            table[index] = newNode;
        }
        else
        {
            Node *currentNode = table[index];
            while (currentNode->next != nullptr)
            {
                currentNode = currentNode->next;
            }
            currentNode->next = newNode;
        }
    }

    int search(const string &key)
    {
        int index = hash(key);

        if (table[index] != nullptr)
        {
            Node *currentNode = table[index];
            while (currentNode != nullptr)
            {
                if (currentNode->key == key)
                {
                    return currentNode->value;
                }
                currentNode = currentNode->next;
            }
        }

        return -1; // Key not found
    }

    void remove(const string &key)
    {
        int index = hash(key);

        if (table[index] != nullptr)
        {
            Node *currentNode = table[index];
            Node *prevNode = nullptr;
            while (currentNode != nullptr)
            {
                if (currentNode->key == key)
                {
                    if (prevNode == nullptr)
                    {
                        table[index] = currentNode->next;
                    }
                    else
                    {
                        prevNode->next = currentNode->next;
                    }
                    delete currentNode;
                    return;
                }
                prevNode = currentNode;
                currentNode = currentNode->next;
            }
        }
    }
};

class Graph
{
private:
    struct Node
    {
        string city;
        int vertex;
        Node *next;
    };

    Node *adjacencyList[MAX_VERTICES];
    int numVertices;

    unordered_map<int, string> cityNames;
    //vector<City> cityList; // Add this line to store the cities

public:
    unordered_map<string, int> cityIndexMap;
     vector<City> cityList; // Add this line to store the cities
    Graph() : numVertices(0)
    {
        for (int i = 0; i < MAX_VERTICES; i++)
        {
            adjacencyList[i] = nullptr;
        }
    }
    int getNumVertices() const
    {
        return numVertices;
    }

    void addEdge(int source, int destination)
    {
        Node *newNode = new Node;
        newNode->vertex = destination;
        newNode->city = cityNameFromIndex(destination);
        newNode->next = nullptr;

        newNode->next = adjacencyList[source];
        adjacencyList[source] = newNode;
    }

    void printNeighbours(string name)
    {
        int cityIndex = cityIndexMap[name];
        if (cityIndexMap.count(name) == 0)
        {
            cout << "City not found." << endl;
            return;
        }

        cout << "Neighboring cities for " << name << ": ";

        Node *currentNode = adjacencyList[cityIndex];
        while (currentNode != nullptr)
        {
            int neighborIndex = currentNode->vertex;
            string neighborName = cityNameFromIndex(neighborIndex);
            cout << neighborName << " ";

            currentNode = currentNode->next;
        }

        cout << endl;
    }

    void addCitiesFromCSV(const string &filename)
    {
        ifstream file(filename);
        if (!file)
        {
            cout << "Failed to open the CSV file." << endl;
            return;
        }

        string line;
        getline(file, line); // Skip the header line

        while (getline(file, line))
        {
            stringstream ss(line);
            string cityName, latitude, longitude, countryName, population;

            getline(ss, cityName, ',');
            getline(ss, latitude, ',');
            getline(ss, longitude, ',');
            getline(ss, countryName, ',');
            getline(ss, population, ',');

            double lat = stod(latitude);
            double lon = stod(longitude);
            int pop;

            try
            {
                size_t pos;
                pop = stoi(population, &pos);

                if (pos < population.length())
                {
                    throw invalid_argument("Invalid population value.");
                }
            }
            catch (const invalid_argument &e)
            {
                cout << "Invalid population value: " << population << endl;
                continue;
            }

            City city{cityName, lat, lon, countryName, pop};
            // Store the city in the cityList vector
            cityList.push_back(city);

            adjacencyList[numVertices] = nullptr;

            cityIndexMap[city.cityName] = numVertices;
            cityNames[numVertices] = city.cityName;

            for (int i = 0; i < numVertices; i++)
            {
                city.neighbors[i] = 1;
                addEdge(numVertices, i);
            }

            numVertices++;
        }

        file.close();
    }

    string cityNameFromIndex(int index)
    {
        if (cityNames.count(index) > 0)
        {
            return cityNames[index];
        }
        else
        {
            return "";
        }
    }

    double calculateDistance(double lat1, double lon1, double lat2, double lon2)
    {
        const double R = 6371.0;

        double lat1Rad = lat1 * M_PI / 180.0;
        double lon1Rad = lon1 * M_PI / 180.0;
        double lat2Rad = lat2 * M_PI / 180.0;
        double lon2Rad = lon2 * M_PI / 180.0;

        double latDiff = lat2Rad - lat1Rad;
        double lonDiff = lon2Rad - lon1Rad;

        double a = sin(latDiff / 2) * sin(latDiff / 2) +
                   cos(lat1Rad) * cos(lat2Rad) *
                       sin(lonDiff / 2) * sin(lonDiff / 2);
        double c = 2 * atan2(sqrt(a), sqrt(1 - a));
        double distance = R * c;

        return distance;
    }

    void bfsTraversal(int startVertex)
    {
        bool visited[MAX_VERTICES] = {false};
        queue<int> q;

        visited[startVertex] = true;
        q.push(startVertex);

        while (!q.empty())
        {
            int currentVertex = q.front();
            q.pop();

            cout << cityNameFromIndex(currentVertex) << " ";

            Node *currentNode = adjacencyList[currentVertex];
            while (currentNode != nullptr)
            {
                int neighborIndex = currentNode->vertex;
                if (!visited[neighborIndex])
                {
                    visited[neighborIndex] = true;
                    q.push(neighborIndex);
                }
                currentNode = currentNode->next;
            }
        }

        cout << endl;
    }
    void dfsTraversal(int startVertex)
    {
        bool visited[MAX_VERTICES] = {false};
        stack<int> s;

        s.push(startVertex);

        while (!s.empty())
        {
            int currentVertex = s.top();
            s.pop();

            if (!visited[currentVertex])
            {
                cout << cityNameFromIndex(currentVertex) << " ";
                visited[currentVertex] = true;
            }

            Node *currentNode = adjacencyList[currentVertex];
            while (currentNode != nullptr)
            {
                int neighborIndex = currentNode->vertex;
                if (!visited[neighborIndex])
                {
                    s.push(neighborIndex);
                }
                currentNode = currentNode->next;
            }
        }

        cout << endl;
    }

    void dijkstra(int startVertex, int destination)
    {
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        vector<double> distance(numVertices, numeric_limits<double>::infinity());
        vector<int> parent(numVertices, -1);
        vector<bool> visited(numVertices, false);

        pq.push(make_pair(0.0, startVertex));
        distance[startVertex] = 0.0;

        while (!pq.empty())
        {
            int currentVertex = pq.top().second;
            pq.pop();

            if (visited[currentVertex])
            {
                continue;
            }

            visited[currentVertex] = true;

            if (currentVertex == destination)
            {
                break;
            }

            Node *temp = adjacencyList[currentVertex];
            while (temp != nullptr)
            {
                int neighbor = temp->vertex;
                double weight = calculateDistance(cityList[currentVertex].latitude, cityList[currentVertex].longitude,
                                                  cityList[neighbor].latitude, cityList[neighbor].longitude);

                if (!visited[neighbor] && distance[currentVertex] + weight < distance[neighbor])
                {
                    distance[neighbor] = distance[currentVertex] + weight;
                    parent[neighbor] = currentVertex;
                    pq.push(make_pair(distance[neighbor], neighbor));
                }

                temp = temp->next;
            }
        }

        if (parent[destination] == -1)
        {
            cout << "Path from " << cityList[startVertex].cityName << " to " << cityList[destination].cityName << " not found." << endl;
        }
        else
        {
            printShortestPath(startVertex, destination, parent);
        }
    }

    void printShortestPath(int startVertex, int destination, const vector<int> &parent)
    {
        stack<int> pathStack;
        int currentVertex = destination;
        pathStack.push(currentVertex);

        while (currentVertex != startVertex)
        {
            currentVertex = parent[currentVertex];
            pathStack.push(currentVertex);
        }

        cout << "Shortest path from " << cityList[startVertex].cityName << " to " << cityList[destination].cityName << ": ";
        while (!pathStack.empty())
        {
            cout << cityList[pathStack.top()].cityName;
            pathStack.pop();
            if (!pathStack.empty())
            {
                cout << " -> ";
            }
        }
        cout << endl;
    }
    void bubbleSort(City *cityList ,int size)
    {
        for (int i = 0; i < size - 1; i++)
        {
            for (int j = 0; j < size - i - 1; j++)
            {
                if (cityList[j].population > cityList[j + 1].population)
                {
                    swap(cityList[j], cityList[j + 1]);
                }
            }
        }
    }


void selectionSort(City *cityList, int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        int minIndex = i;
        for (int j = i + 1; j < size; j++)
        {
            if (cityList[j].population < cityList[minIndex].population)
            {
                minIndex = j;
            }
        }
        swap(cityList[i], cityList[minIndex]);
    }
}



void insertionSort(City *cityList, int size)
{
    for (int i = 1; i < size; i++)
    {
        City key = cityList[i];
        int j = i - 1;
        while (j >= 0 && cityList[j].population > key.population)
        {
            cityList[j + 1] = cityList[j];
            j--;
        }
        cityList[j + 1] = key;
    }
}

void merge(City *cityList, int left, int middle, int right)
{
    int n1 = middle - left + 1;
    int n2 = right - middle;

    City *leftTemp = new City[n1];
    City *rightTemp = new City[n2];

    for (int i = 0; i < n1; i++)
    {
        leftTemp[i] = cityList[left + i];
    }
    for (int i = 0; i < n2; i++)
    {
        rightTemp[i] = cityList[middle + 1 + i];
    }

    int i = 0;
    int j = 0;
    int k = left;

    while (i < n1 && j < n2)
    {
        if (leftTemp[i].population <= rightTemp[j].population)
        {
            cityList[k] = leftTemp[i];
            i++;
        }
        else
        {
            cityList[k] = rightTemp[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        cityList[k] = leftTemp[i];
        i++;
        k++;
    }

    while (j < n2)
    {
        cityList[k] = rightTemp[j];
        j++;
        k++;
    }

    delete[] leftTemp;
    delete[] rightTemp;
}
void mst(const vector<City>& cities, int numVertices)
{
    vector<bool> inMST(numVertices, false);  // Track if a vertex is included in MST
    vector<int> parent(numVertices, -1);     // Track the parent of each vertex in MST
    vector<double> key(numVertices, numeric_limits<double>::infinity());  // Key values to determine the minimum weight

    // Priority queue to store the vertices and their corresponding key values
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    int startVertex = 0;  // Start the MST from vertex 0
    key[startVertex] = 0.0;  // Set the key value of the start vertex to 0
    pq.push(make_pair(0.0, startVertex));  // Push the start vertex into the priority queue

    while (!pq.empty())
    {
        int u = pq.top().second;  // Get the vertex with the minimum key value
        pq.pop();

        inMST[u] = true;  // Mark the vertex as included in MST

        // Iterate through the neighbors of the current vertex
        for (int v : cities[u].neighbors)
        {
            double weight = calculateDistance(cities[u].latitude ,cities[u].longitude , cities[v].latitude ,cities[v].longitude);  // Calculate the weight between the current vertex and its neighbor

            // If the neighbor is not yet included in MST and the weight is smaller than its current key value
            if (!inMST[v] && weight < key[v])
            {
                parent[v] = u;  // Update the parent of the neighbor
                key[v] = weight;  // Update the key value of the neighbor
                pq.push(make_pair(key[v], v));  // Push the neighbor into the priority queue
            }
        }
    }

    // Print the minimum spanning tree
    cout << "Minimum Spanning Tree (MST):" << endl;
    for (int i = 1; i < numVertices; ++i)
    {
        cout << cities[parent[i]].cityName << " - " << cities[i].cityName << endl;
    }
}
};
// Linear search
int linearSearch(const vector<City>& cityList, const string& cityName)
{
    for (int i = 0; i < cityList.size(); i++)
    {
        if (cityList[i].cityName == cityName)
        {
            return i;  // Found at index i
        }
    }
    return -1;  // Not found
}

// Binary search (requires sorted cityList)
int binarySearch(const std::string& cityName, const City cityList[], int size) {
    int low = 0;
    int high = size - 1;

    while (low <= high) {
        int mid = low + (high - low) / 2;
        const std::string& currentCity = cityList[mid].cityName;

        if (currentCity == cityName) {
            return mid;
        } else if (currentCity < cityName) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return -1; // City not found
}
// Hash search using the HashTable class
int hashSearch(HashTable& hashTable, const unordered_map<std::string, int>& cityIndexMap, const string& cityName)
{
    return hashTable.search(cityName);
}



int main()
{
    Graph graph;
    int choice;
    graph.addCitiesFromCSV("C:\\Users\\anwar\\OneDrive\\Desktop\\Dsa Terminal prj\\output\\worldcities.csv");
    cout << "Enter the city name : ";
    string name1, name2;

    cin >> name1;
    cout << "Enter the 2nd city name : ";
    cin >> name2;
    cout << "BFS Traversal of the graph from " << name1 << endl;
    graph.bfsTraversal(graph.cityIndexMap[name1]);
    cout << " " << endl;
    cout << "DFS Traversal of the graph from " << name1 << endl;
    graph.dfsTraversal(graph.cityIndexMap[name1]);
    cout << " " << endl;
    graph.dijkstra(graph.cityIndexMap[name1], graph.cityIndexMap[name2]);
    cout << " " << endl;

    int size = graph.getNumVertices();
    City *cityList = new City[size];

    for (int i = 0; i < size; i++)
    {
        cityList[i] = graph.cityList[i];
    }

    // Sort the cityList using bubble sort
    graph.bubbleSort(cityList, size);
    cout << "Sorted cityList using bubble sort:" << endl;
    for (int i = 0; i < size; i++)
    {
        cout << cityList[i].cityName << ": " << cityList[i].population << endl;
    }
    cout << endl;

    // Sort the cityList using selection sort
    for (int i = 0; i < size; i++)
    {
        cityList[i] = graph.cityList[i];
    }
    graph.selectionSort(cityList, size);
    cout << "Sorted cityList using selection sort:" << endl;
    for (int i = 0; i < size; i++)
    {
        cout << cityList[i].cityName << ": " << cityList[i].population << endl;
    }
    cout << endl;

    // Sort the cityList using insertion sort
    for (int i = 0; i < size; i++)
    {
        cityList[i] = graph.cityList[i];
    }
    graph.insertionSort(cityList, size);
    cout << "Sorted cityList using insertion sort:" << endl;
    for (int i = 0; i < size; i++)
    {
        cout << cityList[i].cityName << ": " << cityList[i].population << endl;
    }
    cout << endl;
    string cityName;
    cout << "Enter the city name: ";
    cin >> cityName;
    
    // Linear search
    int linearIndex = linearSearch(graph.cityList, cityName);
    if (linearIndex != -1)
    {
        cout << "Linear search: City found at index " << linearIndex << endl;
    }
    else
    {
        cout << "Linear search: City not found" << endl;
    }
    
    // Binary search (assuming cityList is sorted)
    //int size = sizeof(cityList) / sizeof(cityList[0]);

    int index = binarySearch(cityName, cityList, size);
    if (index != -1) {
        std::cout << "City: " << cityName << ", Index: " << index << std::endl;
    } else {
        std::cout << "City not found." << std::endl;
    }
    
    // Hash search using the HashTable
    // Assuming graph.cityIndexMap is of type std::unordered_map<std::string, int>
    HashTable hashTable;
    int hashIndex = hashSearch(hashTable, graph.cityIndexMap, cityName);

    if (hashIndex != -1)
    {
        cout << "Hash search: City found at index " << hashIndex << endl;
    }
    else
    {
        cout << "Hash search: City not found" << endl;
    }
    graph.mst(graph.cityList , size);

    delete[] cityList;


    return 0;
}
