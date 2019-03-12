using UnityEngine;
using System.Collections;

public class Time : MonoBehaviour
{
    public static Time instance;
    int index = 0;
	// Use this for initialization

    #region Singleton
    void Awake()
    {
        if (instance != null)
        {
            Debug.LogWarning("Created more than one instance of Time!");
            return;
        }
        instance = this;
    }
    #endregion

	// Update is called once per frame
	void FixedUpdate()
	{
        index++;
	}

    public int GetIndex() {
        return index;
    }
}
