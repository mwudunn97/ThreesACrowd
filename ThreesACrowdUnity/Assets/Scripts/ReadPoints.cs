using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ReadPoints : MonoBehaviour {

    public TextAsset pointFile;
    public Transform personPrefab;
    int num_people;
    int num_points;
    Vector2[,] points;
	// Use this for initialization
	void Start()
	{
        OpenPoints();
        for (int i = 0; i < num_people; i++) {
            Transform person;
            person = Instantiate(personPrefab, new Vector3(i * 2.0F, 0, 0), Quaternion.identity);
            PersonScript ps = person.gameObject.GetComponent<PersonScript>();

            Vector2[] person_points = new Vector2[num_points];
            for (int j = 0; j < num_points; j++) {
                person_points[j] = points[i, j];
            }
            ps.points = person_points;
        }

	}
    private void OpenPoints () {
        string pointText = pointFile.text;
        string[] fLines = pointText.Split('\n');
        string[] dim = fLines[0].Split(' ');
        num_people = int.Parse(dim[0]);
        num_points = int.Parse(dim[1]);
        points = new Vector2[num_people, num_points];

        for (int i = 1; i < fLines.Length; i++) {
            if (fLines[i].Length < 3) {
                break;
            }
            string[] pos = fLines[i].Split(' ');
            int index = i - 1;
            int person_index = index % num_people;
            int point_index = index / num_people;

      
            points[person_index, point_index].x = float.Parse(pos[0]);
            points[person_index, point_index].y = float.Parse(pos[1]);
        }

	}
	
	// Update is called once per frame
	void Update () {
		
	}
}
