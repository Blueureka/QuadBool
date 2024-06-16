#pragma once

class AABB {
public:
    double minX, minY, minZ;
    double maxX, maxY, maxZ;
    double tolerance;

    AABB() : minX(0), minY(0), minZ(0), maxX(0), maxY(0), maxZ(0), tolerance(0.0) {}
    AABB(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, double tolerance = 0.0)
        : minX(minX), minY(minY), minZ(minZ), maxX(maxX), maxY(maxY), maxZ(maxZ), tolerance(tolerance) {}

    // Check if a point is inside the AABB with tolerance
    bool contains(double x, double y, double z) const {
        return (x >= minX - tolerance && x <= maxX + tolerance &&
                y >= minY - tolerance && y <= maxY + tolerance &&
                z >= minZ - tolerance && z <= maxZ + tolerance);
    }

    // Check if this AABB contains another AABB with tolerance
    bool contains(const AABB& other) const {
        return contains(other.minX, other.minY, other.minZ) &&
               contains(other.maxX, other.maxY, other.maxZ);
    }
    
    // Check if two AABBs intersect with tolerance
    bool intersects(const AABB& other) const {
        return (minX <= other.maxX + tolerance && maxX >= other.minX - tolerance &&
                minY <= other.maxY + tolerance && maxY >= other.minY - tolerance &&
                minZ <= other.maxZ + tolerance && maxZ >= other.minZ - tolerance);
    }

    // Calculate the volume of the AABB
    double volume() const {
        return (maxX - minX) * (maxY - minY) * (maxZ - minZ);
    }

    // Calculate the center of the AABB
    void center(double& cx, double& cy, double& cz) const {
        cx = (minX + maxX) / 2;
        cy = (minY + maxY) / 2;
        cz = (minZ + maxZ) / 2;
    }
};

