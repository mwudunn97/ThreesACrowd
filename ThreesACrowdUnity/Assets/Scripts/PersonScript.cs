using UnityEngine;
using System.Collections;

public class PersonScript : MonoBehaviour
{
    public Vector2[] points;

	// Use this for initialization
	void Start()
	{

	}

	// Update is called once per frame
	void Update()
	{
        if (points.Length == 0) {
            return;
        } 
        int index = Time.instance.GetIndex() % points.Length;
        Vector2 point = points[index];
        Vector3 nextLoc = new Vector3(point.x, 0.0f, point.y);
        GameObject plane = GameObject.FindGameObjectWithTag("Plane");
        Bounds bounds = plane.GetComponent<Renderer>().bounds;
        Vector3 nextPos = bounds.min + Vector3.Scale(nextLoc, bounds.size);
        nextPos += new Vector3(0.0f, 0.5f, 0.0f);
        Debug.Log(bounds.size);
        this.gameObject.transform.position = nextPos;
	}
}
