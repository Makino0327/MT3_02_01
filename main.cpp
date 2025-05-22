#include <Novice.h>
#include <cstdint>
#include <imgui.h>
#define _USE_MATH_DEFINES
#include <cmath>


const char kWindowTitle[] = "LC1C_24_マキノハルト_タイトル";

struct Matrix4x4
{
	float m[4][4];
};
struct Vector3
{
	float x, y, z;
};
struct Vector4
{
	float x, y, z, w;
};

struct Line
{
	Vector3 origin;
	Vector3 diff;
};

struct Ray
{
	Vector3 origin;
	Vector3 diff;
};

struct Segment
{
	Vector3 origin;
	Vector3 diff;
};
struct Sphere
{
	Vector3 center;
	float radius;
};
static const int kWindowWidth = 1280;
static const int kWindowHeight = 720;

Vector3 Add(const Vector3& v1, const Vector3& v2) {
	return {
		v1.x + v2.x,
		v1.y + v2.y,
		v1.z + v2.z
	};
}
float Determinant3x3(float matrix[3][3]) {
	return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
		matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
		matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}
float Minor(const Matrix4x4& m, int row, int col) {
	float sub[3][3];
	int sub_i = 0;
	for (int i = 0; i < 4; i++) {
		if (i == row) continue;
		int sub_j = 0;
		for (int j = 0; j < 4; j++) {
			if (j == col) continue;
			sub[sub_i][sub_j] = m.m[i][j];
			sub_j++;
		}
		sub_i++;
	}

	// 3x3行列の行列式を計算
	return Determinant3x3(sub);
}
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return result;
}
Vector3 TransformWithW(const Vector3& v, const Matrix4x4& m) {
	float x = v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + m.m[3][0];
	float y = v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + m.m[3][1];
	float z = v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + m.m[3][2];
	float w = v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + m.m[3][3];

	if (w != 0.0f) {
		x /= w;
		y /= w;
		z /= w;
	}
	return { x, y, z };
}
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	float x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
	float y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
	float z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

	if (w != 0.0f) {
		x /= w;
		y /= w;
		z /= w;
	}

	return { x, y, z };
}
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 mat = {};

	float halfWidth = width * 0.5f;
	float halfHeight = height * 0.5f;
	float depthRange = maxDepth - minDepth;

	mat.m[0][0] = halfWidth;
	mat.m[1][1] = -halfHeight; // Y軸が上から下へ向かう場合
	mat.m[2][2] = depthRange;
	mat.m[3][0] = left + halfWidth;
	mat.m[3][1] = top + halfHeight;
	mat.m[3][2] = minDepth;
	mat.m[3][3] = 1.0f;

	return mat;
}
// x軸回転行列
Matrix4x4 MakeRoteXMatrix(float radian) {
	Matrix4x4 result = {
		1,0,0,0,
		0,cosf(radian),sinf(radian),0,
		0,-sinf(radian),cosf(radian),0,
		0,0,0,1,
	};
	return result;
}
// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result = {
		cosf(radian),0,-sinf(radian),0,
		0,1,0,0,
		sinf(radian),0,cosf(radian),0,
		0,0,0,1,
	};
	return result;
}
// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result = {
		cosf(radian),sinf(radian),0,0,
		-sinf(radian),cosf(radian),0,0,
		0,0,1,0,
		0,0,0,1,
	};
	return result;
}
// 平行移動行列
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = {
	  1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1, 0,
	 translate.x, translate.y, translate.z, 1
	};
	return result;
}
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 scaleMatrix = {
		scale.x, 0, 0, 0,
		0, scale.y, 0, 0,
		0, 0, scale.z, 0,
		0, 0, 0, 1
	};
	Matrix4x4 rotateXMatrix = MakeRoteXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateMatrix = Multiply(Multiply(rotateXMatrix, rotateYMatrix), rotateZMatrix);
	Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);
	return Multiply(Multiply(scaleMatrix, rotateMatrix), translateMatrix);
}
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};

	float f = 1.0f / tanf(fovY / 2.0f);

	result.m[0][0] = f / aspectRatio;
	result.m[1][1] = f;
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1.0f;
	result.m[3][2] = -nearClip * farClip / (farClip - nearClip);
	result.m[3][3] = 0.0f;

	return result;
}
Vector3 Cross(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y = v1.z * v2.x - v1.x * v2.z;
	result.z = v1.x * v2.y - v1.y * v2.x;
	return result;
}
Matrix4x4 MakeLookAtMatrix(const Vector3& eye, const Vector3& target, const Vector3& up) {
	Vector3 zAxis = {
	eye.x - target.x,
	eye.y - target.y,
	eye.z - target.z
	};

	// 正規化
	float lengthZ = sqrtf(zAxis.x * zAxis.x + zAxis.y * zAxis.y + zAxis.z * zAxis.z);
	zAxis = { zAxis.x / lengthZ, zAxis.y / lengthZ, zAxis.z / lengthZ };

	// X軸 = Up × Z
	Vector3 xAxis = Cross(up, zAxis);
	float lengthX = sqrtf(xAxis.x * xAxis.x + xAxis.y * xAxis.y + xAxis.z * xAxis.z);
	xAxis = { xAxis.x / lengthX, xAxis.y / lengthX, xAxis.z / lengthX };

	// Y軸 = Z × X
	Vector3 yAxis = Cross(zAxis, xAxis);

	Matrix4x4 view{};
	view.m[0][0] = xAxis.x;
	view.m[1][0] = xAxis.y;
	view.m[2][0] = xAxis.z;
	view.m[3][0] = -(xAxis.x * eye.x + xAxis.y * eye.y + xAxis.z * eye.z);

	view.m[0][1] = yAxis.x;
	view.m[1][1] = yAxis.y;
	view.m[2][1] = yAxis.z;
	view.m[3][1] = -(yAxis.x * eye.x + yAxis.y * eye.y + yAxis.z * eye.z);

	view.m[0][2] = zAxis.x;
	view.m[1][2] = zAxis.y;
	view.m[2][2] = zAxis.z;
	view.m[3][2] = -(zAxis.x * eye.x + zAxis.y * eye.y + zAxis.z * eye.z);

	view.m[0][3] = 0.0f;
	view.m[1][3] = 0.0f;
	view.m[2][3] = 0.0f;
	view.m[3][3] = 1.0f;

	return view;
}
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);

	// X方向の線
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		float x = -kGridHalfWidth + kGridEvery * xIndex;
		Vector3 worldStart = { x, 0.0f, -kGridHalfWidth };
		Vector3 worldEnd = { x, 0.0f,  kGridHalfWidth };

		Vector3 clipStart = TransformWithW(worldStart, viewProjectionMatrix);
		Vector3 clipEnd = TransformWithW(worldEnd, viewProjectionMatrix);

		Vector3 screenStart = Transform(clipStart, viewportMatrix);
		Vector3 screenEnd = Transform(clipEnd, viewportMatrix);

		uint32_t color = (fabsf(x) < 0.001f) ? 0x000000FF : 0xAAAAAAFF;
		Novice::DrawLine(int(screenStart.x + 0.5f), int(screenStart.y + 0.5f),
			int(screenEnd.x + 0.5f), int(screenEnd.y + 0.5f), color);
	}

	// Z方向の線
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		float z = -kGridHalfWidth + kGridEvery * zIndex;
		Vector3 worldStart = { -kGridHalfWidth, 0.0f, z };
		Vector3 worldEnd = { kGridHalfWidth, 0.0f, z };

		Vector3 clipStart = TransformWithW(worldStart, viewProjectionMatrix);
		Vector3 clipEnd = TransformWithW(worldEnd, viewProjectionMatrix);

		Vector3 screenStart = Transform(clipStart, viewportMatrix);
		Vector3 screenEnd = Transform(clipEnd, viewportMatrix);

		uint32_t color = (fabsf(z) < 0.001f) ? 0x000000FF : 0xAAAAAAFF;
		Novice::DrawLine(int(screenStart.x + 0.5f), int(screenStart.y + 0.5f),
			int(screenEnd.x + 0.5f), int(screenEnd.y + 0.5f), color);
	}
}
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	return {
		v1.x - v2.x,
		v1.y - v2.y,
		v1.z - v2.z
	};
}
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 16;
	const float kLonEvery = 2.0f * float(M_PI) / float(kSubdivision);
	const float kLatEvery = float(M_PI) / float(kSubdivision);

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -float(M_PI) / 2.0f + kLatEvery * (latIndex + 0.5f);

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;

			// ワールド座標での点 a, b, c を計算
			Vector3 a = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon)
			};

			Vector3 b = {
				sphere.center.x + sphere.radius * cosf(lat + kLatEvery) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
				sphere.center.z + sphere.radius * cosf(lat + kLatEvery) * sinf(lon)
			};

			Vector3 c = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon + kLonEvery),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon + kLonEvery)
			};

			// ビュープロジェクション → ビューポート変換まで
			Vector3 screenA = Transform(TransformWithW(a, viewProjectionMatrix), viewportMatrix);
			Vector3 screenB = Transform(TransformWithW(b, viewProjectionMatrix), viewportMatrix);
			Vector3 screenC = Transform(TransformWithW(c, viewProjectionMatrix), viewportMatrix);

			// 線を描画
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenC.x), int(screenC.y), color);
		}
	}
}
Matrix4x4 MakeRotationMatrix(const Vector3& rotate) {
	float cosX = cosf(rotate.x), sinX = sinf(rotate.x);
	float cosY = cosf(rotate.y), sinY = sinf(rotate.y);
	float cosZ = cosf(rotate.z), sinZ = sinf(rotate.z);

	Matrix4x4 m{};

	m.m[0][0] = cosY * cosZ;
	m.m[0][1] = cosY * sinZ;
	m.m[0][2] = -sinY;
	m.m[0][3] = 0.0f;

	m.m[1][0] = sinX * sinY * cosZ - cosX * sinZ;
	m.m[1][1] = sinX * sinY * sinZ + cosX * cosZ;
	m.m[1][2] = sinX * cosY;
	m.m[1][3] = 0.0f;

	m.m[2][0] = cosX * sinY * cosZ + sinX * sinZ;
	m.m[2][1] = cosX * sinY * sinZ - sinX * cosZ;
	m.m[2][2] = cosX * cosY;
	m.m[2][3] = 0.0f;

	m.m[3][0] = 0.0f;
	m.m[3][1] = 0.0f;
	m.m[3][2] = 0.0f;
	m.m[3][3] = 1.0f;

	return m;
}
Matrix4x4 Inverse(const Matrix4x4& m) {
	Matrix4x4 result = {};

	// 4x4行列の行列式を計算
	float det = 0.0f;
	for (int col = 0; col < 4; col++) {
		int sign = (col % 2 == 0) ? 1 : -1;
		det += sign * m.m[0][col] * Minor(m, 0, col);
	}

	// 行列式が0の場合は逆行列が存在しない
	if (det == 0.0f) {
		return result;
	}

	float invDet = 1.0f / det;

	// 各要素の計算
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			int sign = ((i + j) % 2 == 0) ? 1 : -1;
			result.m[j][i] = sign * Minor(m, i, j) * invDet;
		}
	}

	return result;
}
Vector3 Project(const Vector3& v1, const Vector3& v2) {
	// v2の大きさを求める
	float length = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);
	// 正規化ベクトルを求める
	Vector3 unitV2 = { v2.x / length, v2.y / length, v2.z / length };
	// 内積を求める
	float dotProduct = v1.x * unitV2.x + v1.y * unitV2.y + v1.z * unitV2.z;
	// 投影ベクトルを求める
	Vector3 projection = { unitV2.x * dotProduct, unitV2.y * dotProduct, unitV2.z * dotProduct };
	return projection;
}
Vector3 ClosestPoint(const Vector3& point, const Segment& segment) {
	// 線分の長さを求める
	float length = sqrt(segment.diff.x * segment.diff.x + segment.diff.y * segment.diff.y + segment.diff.z * segment.diff.z);
	// 線分の単位ベクトルを求める
	Vector3 unitDiff = { segment.diff.x / length, segment.diff.y / length, segment.diff.z / length };
	// 始点からのベクトルを求める
	Vector3 originToPoint = { point.x - segment.origin.x, point.y - segment.origin.y, point.z - segment.origin.z };
	// 投影ベクトルを求める
	Vector3 projection = Project(originToPoint, unitDiff);
	// 最近接点を求める
	Vector3 closestPoint = { segment.origin.x + projection.x, segment.origin.y + projection.y, segment.origin.z + projection.z };
	return closestPoint;
}
bool IsCollision(const Sphere& sphere, const Sphere& sphere2) {
	// 球の中心間の距離を求める
	Vector3 diff = Subtract(sphere.center, sphere2.center);
	float distance = sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
	// 半径の和を求める
	float radiusSum = sphere.radius + sphere2.radius;
	// 衝突判定
	return distance <= radiusSum;
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	// カメラ
	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	// 球
	Sphere Sphere1 = { { 0.0f, 0.0f, 0.0f }, 0.4f };
	Sphere Sphere2 = { { 0.0f, 0.0f, 1.0f }, 0.4f };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		// グリッド
		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, { 0.0f,0.0f,0.0f });
		Matrix4x4 cameraRotationMatrix = Multiply(Multiply(MakeRoteXMatrix(cameraRotate.x), MakeRotateYMatrix(cameraRotate.y)), MakeRotateZMatrix(cameraRotate.z));
		Matrix4x4 cameraMatrix = Multiply(cameraRotationMatrix, MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, cameraTranslate));
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("Sphere1", &Sphere1.center.x, 0.01f);
		ImGui::DragFloat("Sphere1Radius", &Sphere1.radius, 0.01f);
		ImGui::DragFloat3("Sphere2", &Sphere2.center.x, 0.01f);
		ImGui::DragFloat("Sphere2Radius", &Sphere2.radius, 0.01f);

		ImGui::End();


		DrawGrid(worldViewProjectionMatrix, viewportMatrix);
		if (IsCollision(Sphere1, Sphere2))
		{
			DrawSphere(Sphere1, worldViewProjectionMatrix, viewportMatrix, RED);
		}
		else
		{
			DrawSphere(Sphere1, worldViewProjectionMatrix, viewportMatrix, WHITE);
		}
		DrawSphere(Sphere2, worldViewProjectionMatrix, viewportMatrix, WHITE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
