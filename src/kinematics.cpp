#include "kinematics.h"

#include <algorithm>

#include "utils.h"
void forwardKinematics(const Posture& posture, Bone* bone) {
  // TODO (FK)
  // Same as HW2, but have some minor change
  // Hint:
  //   1. If you don't use `axis` in this function, you can copy-paste your code
  // Note:
  //   1. bone.axis becomes quaternion instead of vector3f
  while (bone != nullptr) {
    if (bone->parent == nullptr) {
      bone->startPosition = bone->endPosition = posture.translations[0];
      bone->rotation = posture.rotations[0];
    } else {
      bone->rotation = bone->parent->rotation * bone->rotationParentCurrent * posture.rotations[bone->idx];
      bone->startPosition = bone->parent->endPosition + posture.translations[bone->idx];
      bone->endPosition = bone->startPosition + bone->rotation * bone->direction * bone->length;
    }
    forwardKinematics(posture, bone->child);
    bone = bone->sibling;
  }
}

Eigen::VectorXf leastSquareSolver(const Eigen::Matrix3Xf& jacobian, const Eigen::Vector3f& target) {
  // TODO (find x which min(| jacobian * x - target |))
  // Hint:
  //   1. Linear algebra - least squares solution
  //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
  // Note:
  //   1. SVD or other pseudo-inverse method is useful
  //   2. Some of them have some limitation, if you use that method you should check it.
  return jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(target);
}

void inverseKinematics(const Eigen::Vector3f& target, Bone* start, Bone* end, Posture& posture) {
  constexpr int maxIterations = 10000;
  constexpr float epsilon = 1E-3f;
  constexpr float step = 0.1f;
  // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is root.
  Bone* root = start - start->idx;
  std::vector<Bone*> boneList;
  // TODO
  // Hint:
  //   1. Traverse from end to start is easier than start to end (since there is only 1 parent)
  //   2. If start bone is not reachable from end. Go to root first.
  // Note:
  //   1. Both start and end should be in the list
  Bone* current = end;
  while (current != nullptr) {
    boneList.emplace_back(current);
    if (current == start) break;
    current = current->parent;
  }
  if (current == nullptr) {
    current = start;
    while (current != nullptr) {
      if (current->idx) boneList.emplace_back(current);
      current = current->parent;
    }
  }
  size_t boneNum = boneList.size();
  Eigen::Matrix3Xf jacobian(3, 3 * boneNum);
  jacobian.setZero();

  for (int i = 0; i < maxIterations; ++i) {
    forwardKinematics(posture, root);
    Eigen::Vector3f desiredVector = target - end->endPosition;
    if (desiredVector.norm() < epsilon) break;
    // TODO (compute jacobian)
    //   1. Compute arm vectors
    //   2. Compute jacobian columns
    // Hint:
    //   1. You should not put rotation in jacobian if it doesn't have that DoF.
    //   2. jacobian.col(/* some column index */) = /* jacobian column */
    for (size_t j = 0; j < boneNum; j++) {
      const auto& bone = *boneList[j];
      Eigen::Vector3f armVector = end->endPosition - bone.startPosition;
      Eigen::Matrix3f rotation = bone.rotation.toRotationMatrix();
      if (bone.dofrx) jacobian.col(3 * j) = rotation.col(0).normalized().cross(armVector);
      if (bone.dofry) jacobian.col(3 * j + 1) = rotation.col(1).normalized().cross(armVector);
      if (bone.dofrz) jacobian.col(3 * j + 2) = rotation.col(2).normalized().cross(armVector);
    }

    Eigen::VectorXf deltaTheta = step * leastSquareSolver(jacobian, desiredVector);
    // TODO (update rotation)
    //   1. Update posture's eulerAngle using deltaTheta
    // Hint:
    //   1. Use posture.eulerAngle to get posture's eulerAngle
    //   2. All angles are in radians.
    //   3. You can ignore rotation limit of the bone.
    // Bonus:
    //   1. You cannot ignore rotation limit of the bone.

    for (size_t j = 0; j < boneNum; j++) {
      const auto& bone = *boneList[j];
      posture.eulerAngle[bone.idx] += deltaTheta.segment<3>(3 * j);
      posture.eulerAngle[bone.idx][0] = std::clamp(posture.eulerAngle[bone.idx][0], bone.rxmin, bone.rxmax);
      posture.eulerAngle[bone.idx][1] = std::clamp(posture.eulerAngle[bone.idx][1], bone.rymin, bone.rymax);
      posture.eulerAngle[bone.idx][2] = std::clamp(posture.eulerAngle[bone.idx][2], bone.rzmin, bone.rzmax);

      posture.rotations[bone.idx] = Eigen::AngleAxisf(posture.eulerAngle[bone.idx][2], Eigen::Vector3f::UnitZ()) *
                                    Eigen::AngleAxisf(posture.eulerAngle[bone.idx][1], Eigen::Vector3f::UnitY()) *
                                    Eigen::AngleAxisf(posture.eulerAngle[bone.idx][0], Eigen::Vector3f::UnitX());
    }
  }
}