using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.005f; // default 0.015
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;                 // for collision

	Mesh mesh;
	Vector3[] vertices;

	float bounce = 0.5f;
	float friction = 0.2f;

	// Use this for initialization
	void Start () 
	{		
		mesh = GetComponent<MeshFilter>().mesh;
		vertices = mesh.vertices;

		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
// 转动惯量的公式为 I_ref = 求和Massi(riT*ri*I - ri*riT)
// riT*ri*I相当于单位矩阵乘上一个ri的平方，这里diag就是向量长度的平方                                  
// 而第二项，举个例子 riT = (1, 0, -1);	
// vertices[i][0]*vertices[i][0] 即为 1 * 1，这里的vertices[i]仅仅是获取ri

			mass += m;
			float diag = m * vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

	// In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		int collisionVerticesNum = 0;
		Vector3 collisionPositions = Vector3.zero;
		float distance = 0;
		Vector3 rotation = transform.eulerAngles;
		Vector3 centerPosition = transform.position;
		foreach (Vector3 v in vertices)
        {
			// 点到平面的距离   
			Vector3 OP = v - P;
			Vector3 PN = P + N;
			float projectionLength = Vector3.Dot(OP, PN);
			// 如果距离大于0就计算下一个点
			if (projectionLength > 0)   
				continue;
			//Vector3 velocityi = this.v + Vector3.Cross(this.w, ) Vector3.Dot(this.w, rotation) * (v - centerPosition);
			//float velocityDotN = Vector3.Dot(velocityi, N);
			// 如果vi dot N 小于0那证明速度已经更新过，就计算下一个点
			//if (velocityDotN < 0)
			//	continue;
			collisionVerticesNum++;
			float dis = Mathf.Sqrt(projectionLength * projectionLength + OP.sqrMagnitude);
			collisionPositions += v;
			distance += dis;	
        }
		// 计算vi new
		if (collisionVerticesNum == 0)
			return;
		collisionPositions /= (float)collisionVerticesNum;
		distance /= (float)collisionVerticesNum;
		Vector3 averageVelocityi = this.v + Vector3.Dot(this.w, rotation) * (collisionPositions - centerPosition);
		Vector3 velocityNi = Vector3.Dot(averageVelocityi, N) * N;
		Vector3 velocityTi = averageVelocityi - velocityNi;
		float parameter_a = Mathf.Max(1 - friction * (1 + bounce) * velocityNi.magnitude / velocityTi.magnitude, 0);
		Vector3 velocityNiNew = -bounce * velocityNi;
		Vector3 velocityTiNew = parameter_a * velocityTi;
		Vector3 velocityiNew = velocityNiNew + velocityTiNew;
		// 计算 impulse j
		//float parameter_K = 1 / mass - 

	}

	// Update is called once per frame
	void Update () 
	{
		//Game Control
		if(Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);

			restitution = 0.5f;
			// 恢复系数，0-1，为1时是完全弹性碰撞，为0时为完全非弹性碰撞

			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
			launched=true;
		}

		// Part I: Update velocities
		v += dt * Physics.gravity;
		Debug.Log(v);
		// Part II: Collision Impulse
		Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		if (!launched)
			return;
		//Update linear status
		Vector3 x    = transform.position;
		x += dt * v;
		//Update angular status
		Quaternion q = transform.rotation;
		Vector3 rotation = q.eulerAngles;
		rotation += dt * w;
		q.eulerAngles = rotation;
		// TODO: 弄清楚那里是不是叉乘
		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
