#include <bits/stdc++.h>
using namespace std;

class Vertex;
class Edge;
class Face;
vector<Vertex *> total_notches;
vector<pair<Vertex *, Vertex *>> adding_edge;

class Vertex
{
public:
    double x, y;
};

class Edge
{
public:
    Vertex *origin = NULL;
    Edge *twin = NULL;
    // Face *face = NULL;
    Edge *prev = NULL;
    Edge *next = NULL;
    Vertex *dest = NULL;
};

class Face
{
public:
    vector<Vertex *> vertices = {};
};

class DCEL
{
public:
    vector<Vertex *> Vert;
    vector<Edge *> Edges;
    vector<Face *> faces;
    map<Vertex *, pair<Edge *, Vertex *>> mp;
    DCEL(vector<Vertex *> &vertices)
    {
        int numVertices = vertices.size();
        Vert = vertices;
        for (int i = 0; i < numVertices; i++)
        {
            Edge *edge = new Edge();
            edge->origin = vertices[i];
            edge->dest = vertices[(i + 1) % numVertices];
            Edge *twin1 = new Edge();
            twin1->origin = vertices[(i + 1) % numVertices];
            twin1->dest = vertices[i];
            twin1->twin = edge;
            edge->twin = twin1;
            Edges.push_back(edge);
        }
        for (int i = 0; i < numVertices; i++)
        {
            Edges[i]->next = Edges[(i + 1) % numVertices];
            Edges[i]->prev = Edges[(i - 1 + numVertices) % numVertices];
            Edges[i]->twin->next = Edges[(i - 1 + numVertices) % numVertices]->twin;
            Edges[i]->twin->prev = Edges[(i + 1) % numVertices]->twin;
        }
        for (int i = 0; i < numVertices; i++)
        {
            mp[vertices[i]] = {Edges[i], vertices[(i + 1) % numVertices]};
        }
    }
};

vector<Vertex *> calc_notches(DCEL *polygon)
{
    int n = polygon->Vert.size();
    vector<Vertex *> ans;
    for (int i = 0; i < n; i++)
    {
        Vertex *v1 = polygon->Vert[(i - 1 + n) % n];
        Vertex *v2 = polygon->Vert[i];
        Vertex *v3 = polygon->Vert[(i + 1) % n];

        if ((v2->x - v1->x) * (v3->y - v2->y) - (v3->x - v2->x) * (v2->y - v1->y) > 0)
        {
            ans.push_back(v2);
        }
    }
    return ans;
}

bool is_notch(Vertex *v1, Vertex *v2, Vertex *v3)
{
    if ((((v2->x - v1->x) * (v3->y - v2->y)) - ((v3->x - v2->x) * (v2->y - v1->y))) <= 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

DCEL *build_polygon(vector<Vertex *> &vertices)
{
    DCEL *polygon = new DCEL(vertices);
    return polygon;
}
vector<double> calc_rectangle(vector<Vertex *> L)
{
    vector<double> rec(4);
    rec[0] = INT_MIN;
    rec[1] = INT_MAX;
    rec[2] = INT_MIN;
    rec[3] = INT_MAX;
    for (int i = 0; i < L.size(); i++)
    {
        if (L[i]->x > rec[0])
        {
            rec[0] = L[i]->x;
        }
        if (L[i]->x < rec[1])
        {
            rec[1] = L[i]->x;
        }
        if (L[i]->y > rec[2])
        {
            rec[2] = L[i]->y;
        }
        if (L[i]->y < rec[3])
        {
            rec[3] = L[i]->y;
        }
    }
    return rec;
}

bool on_right(Vertex *v1, Vertex *notch, Vertex *curr)
{
    double a = notch->x - v1->x;
    double b = notch->y - v1->y;
    double c = curr->x - v1->x;
    double d = curr->y - v1->y;

    if (a * d - b * c <= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool rectangle(Vertex *v, vector<double> rec)
{
    if (v->x <= rec[0] && v->x >= rec[1] && v->y <= rec[2] && v->y >= rec[3])
    {
        return true;
    }
    else
    {
        return false;
    }
}

void solve(DCEL *&polygon)
{
    vector<Vertex *> L;
    vector<Vertex *> list = polygon->Vert;
    Vertex *init_vertex = list[0];
    int index = 0;
    vector<Vertex *> LPVS = calc_notches(polygon);
    Vertex *v1 = init_vertex;
    Vertex *v2 = polygon->mp[v1].second;
    while (list.size() > 3)
    {
        L.push_back(init_vertex);
        v1 = init_vertex;
        v2 = polygon->mp[v1].second;
        L.push_back(v2);

        // 3.3
        while (!is_notch(L[L.size() - 2], L[L.size() - 1], polygon->mp[L[L.size() - 1]].second) && !is_notch(L[L.size() - 1], polygon->mp[L[L.size() - 1]].second, v1) && !is_notch(polygon->mp[L[L.size() - 1]].second, v1, v2) && L.size() < list.size())
        {
            L.push_back(polygon->mp[L[L.size() - 1]].second);
        }

        // 3.4
        if (list.size() != L.size())
        {
            vector<Vertex *> notches;

            for (int i = 0; i < total_notches.size(); i++)
            {
                int flag = 0;
                for (int j = 0; j < L.size(); j++)
                {
                    if (total_notches[i]->x == L[j]->x && total_notches[i]->y == L[j]->y)
                    {
                        flag--;
                    }
                }
                for (int j = 0; j < list.size(); j++)
                {
                    if (total_notches[i]->x == list[j]->x && total_notches[i]->y == list[j]->y)
                    {
                        flag++;
                        break;
                    }
                }

                if (flag == 1)
                {
                    notches.push_back(total_notches[i]);
                }
            }
            vector<double> rec = calc_rectangle(L);
            for (int i = 0; i < notches.size(); i++)
            {
                if (!rectangle(notches[i], rec))
                {
                    notches.erase(notches.begin() + i, notches.begin() + i + 1);
                    i--;
                }
            }

            for (int i = 0; i < notches.size(); i++)
            {
                // added this
                if (!on_right(v1, notches[i], L[L.size() - 1]))
                {
                    continue;
                }

                for (int j = 1; j < L.size() && L.size() > 2; j++)
                {
                    bool right = on_right(v1, notches[i], L[j]);
                    if (right == true)
                    {
                        L.erase(L.begin() + j, L.begin() + j + 1);
                        j--;
                    }
                }
            }
        }

        // 3.5
        if (L.size() >= 3)
        // if(L[L.size()-1]!=v2)
        {
            // Create Faces
            Edge *new_edge = new Edge();
            new_edge->next = polygon->mp[L[L.size() - 1]].first;
            polygon->mp[L[0]] = {new_edge, L[L.size() - 1]};

            pair<Vertex *, Vertex *> test = {L[0], L[L.size() - 1]};
            adding_edge.push_back(test);
            cout << L[0]->x << " " << L[0]->y << " || " << L[L.size() - 1]->x << " " << L[L.size() - 1]->y << endl;
            // prev, twin, org, dest

            for (int i = 1; i + 1 < L.size(); i++)
            {
                for (int j = 0; j < list.size(); j++)
                {
                    if (L[i]->x == list[j]->x && L[i]->y == list[j]->y)
                    {
                        list.erase(list.begin() + j, list.begin() + j + 1);

                        j--;
                    }
                }
            }
            Face *new_face = new Face();
            for (int i = 0; i < L.size(); i++)
            {
                new_face->vertices.push_back(L[i]);
            }
            polygon->faces.push_back(new_face);

            init_vertex = L[L.size() - 1];
            L.clear();
        }
        else
        {
            init_vertex = L[L.size() - 1];
            L.clear();
        }
    }

    if (list.size() == 3)
    {
        Face *new_face = new Face();
        for (int i = 0; i < list.size(); i++)
        {
            new_face->vertices.push_back(list[i]);
        }
        polygon->faces.push_back(new_face);
    }
}

int main()
{

    freopen("input.txt", "r", stdin);
    int n;
    cin >> n;
    vector<Vertex *> vertices;
    for (int i = 0; i < n; i++)
    {
        Vertex *v = new Vertex();
        double w, u;
        cin >> w >> u;
        v->x = w;
        v->y = u;
        vertices.push_back(v);
    }
    reverse(vertices.begin(), vertices.end());
    DCEL *polygon = build_polygon(vertices);

    vector<Vertex *> t = calc_notches(polygon);
    for (int i = 0; i < t.size(); i++)
    {
        total_notches.push_back(t[i]);
    }

    solve(polygon);

    cout << "FACES" << endl;
    cout << polygon->faces.size() << endl;

    for (int i = 0; i < polygon->faces.size(); i++)
    {
        cout << "FACE"
             << " " << i + 1 << endl;
        for (int j = 0; j < polygon->faces[i]->vertices.size(); j++)
        {
            cout << polygon->faces[i]->vertices[j]->x << " " << polygon->faces[i]->vertices[j]->y << endl;
        }
    }

    std::ofstream outfile("Coordinates.txt");
    if (outfile.is_open())
    {
        for (int i = 0; i < polygon->faces.size(); i++)
        {
            for (int j = 0; j < polygon->faces[i]->vertices.size(); j++)
            {
                outfile << polygon->faces[i]->vertices[j]->x << " " << polygon->faces[i]->vertices[j]->y << " ";
            }
            outfile << "\n";
        }

        outfile.close();
        cout << "Successfully wrote face coordinates to file." << std::endl;
    }
    else
    {
        std::cerr << "Error: Unable to open file." << std::endl;
    }
    return 0;
}