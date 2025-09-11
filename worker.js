// Replaced importScripts with a modern ES6 module import.
import * as THREE from 'https://cdn.skypack.dev/three@0.132.2';

// --- Default Generation Configuration (mirrored from main script) ---
const config = {
    C_C_BOND_LENGTH: 1.42,
    IDEAL_BOND_ANGLE: 2 * Math.PI / 3, // 120 degrees in radians
    BOND_THRESHOLD_FACTOR: 1.2,
    MIN_BOND_LENGTH_FACTOR: 0.8,
};

// This is the main entry point for the worker. It listens for messages from the main thread.
self.onmessage = function(e) {
    const params = e.data;

    // --- Start the multi-step model generation process ---
    
    // 1. Generate the initial lattice
    const packingFactor = 1 - (params.packingDensity / 100);
    const bondThreshold = config.C_C_BOND_LENGTH * config.BOND_THRESHOLD_FACTOR * packingFactor;
    const unitCellDimension = 4;
    const numCells = params.modelSize * unitCellDimension;
    let initialAtoms = generateGrapheneLattice(numCells, numCells, packingFactor, params.rotationAngle);

    // 2. Introduce defects
    let defectedAtoms = introduceStoneWalesDefects(initialAtoms, params.numDefects, bondThreshold);
    
    // 3. Add initial 3D displacement
    let displacedAtoms = addInitialZDisplacement(defectedAtoms, params.zStrength);

    // 4. Determine and refine the bond network to be 3-coordinated
    const atomPositions = displacedAtoms.map(d => d.pos);
    let bondNetwork = determineBondNetwork(atomPositions);
    bondNetwork = refineBondNetwork(atomPositions, bondNetwork);
    
    // 5. Identify the atoms at the edges of the model
    const edgeAtoms = findEdgeAtoms(displacedAtoms);

    // 6. Analyze the ring structure to find defects
    const analysisResult = analyzeRings(bondNetwork, displacedAtoms.length, edgeAtoms);

    // 7. Identify atoms that are part of non-hexagonal rings
    const defectiveAtoms = new Set();
    analysisResult.bondsInNonSixRings.forEach(bondKey => {
        const [u, v] = bondKey.split('-').map(Number);
        defectiveAtoms.add(u);
        defectiveAtoms.add(v);
    });

    // 8. Relax the structure using a physics-based simulation
    let relaxedAtoms = relaxStructure(
        displacedAtoms, 
        bondNetwork, 
        params.relaxIterations, 
        params.lengthStrength, 
        params.angleStrength, 
        params.pristineAngleStrength, 
        params.planarityStrength, 
        defectiveAtoms, 
        config.C_C_BOND_LENGTH, 
        config.IDEAL_BOND_ANGLE
    );

    // --- Send the final, computed data back to the main thread for rendering ---
    self.postMessage({
        type: 'result',
        data: {
            relaxedAtoms: relaxedAtoms,
            bondNetwork: bondNetwork,
            analysisResult: analysisResult,
            edgeAtoms: Array.from(edgeAtoms) // Convert Set to Array for transport
        }
    });
};

// --- All computation functions are now located in the worker ---

function postStatus(message) {
    self.postMessage({ type: 'status', message: message });
}

function generateGrapheneLattice(numCellsX, numCellsY, packingFactor, angleInDegrees) {
    postStatus('Generating perfect graphene lattice...');
    const atomMap = new Map();
    const atom_dist = config.C_C_BOND_LENGTH * packingFactor;
    
    const a1 = new THREE.Vector3(atom_dist * 3 / 2, atom_dist * Math.sqrt(3) / 2, 0);
    const a2 = new THREE.Vector3(atom_dist * 3 / 2, -atom_dist * Math.sqrt(3) / 2, 0);
    
    const b1 = new THREE.Vector3(0, 0, 0);
    const b2 = new THREE.Vector3(atom_dist, 0, 0);

    const baseUnitCellAtoms = [];
    const unitSize = 4;
    const unitCellBox = new THREE.Box3();
    for (let i = 0; i < unitSize; i++) {
        for (let j = 0; j < unitSize; j++) {
            const R = new THREE.Vector3().addScaledVector(a1, i).addScaledVector(a2, j);
            const p1 = new THREE.Vector3().addVectors(R, b1);
            const p2 = new THREE.Vector3().addVectors(R, b2);
            baseUnitCellAtoms.push(p1);
            baseUnitCellAtoms.push(p2);
            unitCellBox.expandByPoint(p1);
            unitCellBox.expandByPoint(p2);
        }
    }
    const unitCellCenter = unitCellBox.getCenter(new THREE.Vector3());

    const rotationAngle = angleInDegrees * (Math.PI / 180);
    const rotationAxis = new THREE.Vector3(0, 0, 1);
    const translationVectorX = new THREE.Vector3().addScaledVector(a1, unitSize);
    const translationVectorY = new THREE.Vector3().addScaledVector(a2, unitSize);
    const numTranslationsX = Math.ceil(numCellsX / unitSize);
    const numTranslationsY = Math.ceil(numCellsY / unitSize);

    for (let m = -numTranslationsX; m <= numTranslationsX; m++) {
        for (let n = -numTranslationsY; n <= numTranslationsY; n++) {
            const totalTranslation = new THREE.Vector3().addScaledVector(translationVectorX, m).addScaledVector(translationVectorY, n);
            const isRotatedTile = (m + n) % 2 !== 0;

            if (isRotatedTile) {
                const tileCenter = new THREE.Vector3().addVectors(totalTranslation, unitCellCenter);
                baseUnitCellAtoms.forEach(atomLocalPos => {
                    const posRelativeToUnitCenter = new THREE.Vector3().subVectors(atomLocalPos, unitCellCenter);
                    posRelativeToUnitCenter.applyAxisAngle(rotationAxis, rotationAngle);
                    const finalPos = new THREE.Vector3().addVectors(posRelativeToUnitCenter, tileCenter);
                    const key = `${finalPos.x.toFixed(4)},${finalPos.y.toFixed(4)},${finalPos.z.toFixed(4)}`;
                    if (!atomMap.has(key)) {
                        atomMap.set(key, { pos: [finalPos.x, finalPos.y, finalPos.z], tile: { m, n } });
                    }
                });
            } else {
                baseUnitCellAtoms.forEach(atomPos => {
                    const newPos = new THREE.Vector3().addVectors(atomPos, totalTranslation);
                    const key = `${newPos.x.toFixed(4)},${newPos.y.toFixed(4)},${newPos.z.toFixed(4)}`;
                    if (!atomMap.has(key)) {
                        atomMap.set(key, { pos: [newPos.x, newPos.y, newPos.z], tile: { m, n } });
                    }
                });
            }
        }
    }
    return Array.from(atomMap.values());
}

function introduceStoneWalesDefects(atomData, defectCount, threshold) {
    postStatus('Introducing defects...');
    if (defectCount === 0) return atomData;
    
    let positions = atomData.map(a => new THREE.Vector3(...a.pos));
    let createdDefects = 0;
    let totalAttempts = 0;
    const maxTotalAttempts = defectCount * 200; 

    while (createdDefects < defectCount && totalAttempts < maxTotalAttempts) {
        totalAttempts++;
        const adj = new Map(positions.map((_, i) => [i, []]));
        const bonds = [];
        for (let i = 0; i < positions.length; i++) {
            for (let j = i + 1; j < positions.length; j++) {
                if (positions[i].distanceTo(positions[j]) < threshold) {
                    adj.get(i).push(j);
                    adj.get(j).push(i);
                    bonds.push([i, j]);
                }
            }
        }
        
        if (bonds.length === 0) break;

        const randomBond = bonds[Math.floor(Math.random() * bonds.length)];
        const [i, j] = randomBond;
        if (adj.get(i).length !== 3 || adj.get(j).length !== 3) continue;

        const p_i = positions[i], p_j = positions[j];
        const midpoint = new THREE.Vector3().addVectors(p_i, p_j).multiplyScalar(0.5);
        const rotationAxis = new THREE.Vector3(0, 0, 1);
        const p_i_new = p_i.clone().sub(midpoint).applyAxisAngle(rotationAxis, Math.PI / 2).add(midpoint);
        const p_j_new = p_j.clone().sub(midpoint).applyAxisAngle(rotationAxis, Math.PI / 2).add(midpoint);

        let clashes = false;
        const minDistance = config.C_C_BOND_LENGTH * config.MIN_BOND_LENGTH_FACTOR;
        for (let k = 0; k < positions.length; k++) {
            if (k === i || k === j) continue;
            if (p_i_new.distanceTo(positions[k]) < minDistance || p_j_new.distanceTo(positions[k]) < minDistance) {
                clashes = true;
                break;
            }
        }

        if (!clashes) {
            positions[i].copy(p_i_new);
            positions[j].copy(p_j_new);
            createdDefects++;
            if (createdDefects % 20 === 0) postStatus(`Introducing defects... (${createdDefects}/${defectCount})`);
        }
    }
    
    return atomData.map((data, i) => ({ pos: [positions[i].x, positions[i].y, positions[i].z], tile: data.tile }));
}

function addInitialZDisplacement(atomData, zStrength) {
    if (zStrength === 0) return atomData;
    postStatus('Applying initial Z-displacement...');
    return atomData.map(data => ({
        pos: [data.pos[0], data.pos[1], (Math.random() - 0.5) * zStrength],
        tile: data.tile
    }));
}

function determineBondNetwork(atoms) {
    postStatus('Determining initial bond network...');
    const bonds = new Set();
    const atomVectors = atoms.map(a => new THREE.Vector3(...a));
    for (let i = 0; i < atomVectors.length; i++) {
        const neighbors = [];
        for (let j = 0; j < atomVectors.length; j++) {
            if (i === j) continue;
            neighbors.push({ index: j, distance: atomVectors[i].distanceTo(atomVectors[j]) });
        }
        neighbors.sort((a, b) => a.distance - b.distance);
        for (const neighbor of neighbors.slice(0, 3)) {
            bonds.add(`${Math.min(i, neighbor.index)}-${Math.max(i, neighbor.index)}`);
        }
    }
    return Array.from(bonds).map(b => b.split('-').map(Number));
}

function refineBondNetwork(atomPositions, initialBonds) {
    postStatus('Refining bond network...');
    const atomVectors = atomPositions.map(p => new THREE.Vector3(...p));
    const bonds = new Set(initialBonds.map(b => `${Math.min(b[0], b[1])}-${Math.max(b[0], b[1])}`));
    const maxIterations = 10;

    const findAndBreakSmallestRingBond = (currentBonds, adj) => {
        // This function is complex and was causing issues.
        // For stability, we will now use a simpler refinement strategy below
        // that focuses on fixing coordination numbers directly.
        // We will break the longest bond in any found 3 or 4 membered rings.
        const bondsToRemove = new Set();

        // Find and break longest bond in 3-membered rings
        for (let i = 0; i < atomVectors.length; i++) {
            const neighbors = adj.get(i);
            if(neighbors.length < 2) continue;
            for(let j_idx=0; j_idx < neighbors.length; j_idx++) {
                for(let k_idx=j_idx+1; k_idx < neighbors.length; k_idx++) {
                    const j = neighbors[j_idx];
                    const k = neighbors[k_idx];
                    if (adj.get(j).includes(k)) { // Found a 3-ring
                         const ringBonds = [
                            {atoms: [i,j], dist: atomVectors[i].distanceTo(atomVectors[j])},
                            {atoms: [j,k], dist: atomVectors[j].distanceTo(atomVectors[k])},
                            {atoms: [k,i], dist: atomVectors[k].distanceTo(atomVectors[i])}
                         ];
                         ringBonds.sort((a,b)=>b.dist - a.dist);
                         const bondToRemove = ringBonds[0].atoms;
                         bondsToRemove.add(`${Math.min(bondToRemove[0], bondToRemove[1])}-${Math.max(bondToRemove[0], bondToRemove[1])}`);
                    }
                }
            }
        }
        // If we found 3-rings to break, do that and return.
        if (bondsToRemove.size > 0) {
            bondsToRemove.forEach(b => currentBonds.delete(b));
            return true;
        }

        // Find and break longest bond in 4-membered rings
        for(let i=0; i<atomVectors.length; i++) {
             const neighborsOfI = adj.get(i);
             if (neighborsOfI.length < 2) continue;
             for (let j_idx = 0; j_idx < neighborsOfI.length; j_idx++) {
                 for(let l_idx = j_idx + 1; l_idx < neighborsOfI.length; l_idx++) {
                    const j = neighborsOfI[j_idx];
                    const l = neighborsOfI[l_idx];
                    const commonNeighbors = adj.get(j).filter(n => n !== i && adj.get(l).includes(n));
                    for(const k of commonNeighbors) {
                         const ringBonds = [
                             {atoms: [i,j], dist: atomVectors[i].distanceTo(atomVectors[j])},
                             {atoms: [j,k], dist: atomVectors[j].distanceTo(atomVectors[k])},
                             {atoms: [k,l], dist: atomVectors[k].distanceTo(atomVectors[l])},
                             {atoms: [l,i], dist: atomVectors[l].distanceTo(atomVectors[i])}
                         ];
                         ringBonds.sort((a,b)=>b.dist - a.dist);
                         const bondToRemove = ringBonds[0].atoms;
                         bondsToRemove.add(`${Math.min(bondToRemove[0], bondToRemove[1])}-${Math.max(bondToRemove[0], bondToRemove[1])}`);
                    }
                 }
             }
        }

        if (bondsToRemove.size > 0) {
            bondsToRemove.forEach(b => currentBonds.delete(b));
            return true;
        }

        return false;
    }

    for (let iter = 0; iter < maxIterations; iter++) {
        let changesMade = false;
        
        const adj = new Map(atomVectors.map((_, i) => [i, []]));
        for(const bondKey of bonds) {
            const [u,v] = bondKey.split('-').map(Number);
            adj.get(u).push(v);
            adj.get(v).push(u);
        }

        if (findAndBreakSmallestRingBond(bonds, adj)) {
            changesMade = true;
            // Rebuild adjacency map after breaking bonds
            for(const [key, val] of adj) val.length = 0;
            for(const bondKey of bonds) {
                const [u,v] = bondKey.split('-').map(Number);
                adj.get(u).push(v);
                adj.get(v).push(u);
            }
        }

        for(let i=0; i<atomVectors.length; i++) {
            if (adj.get(i).length > 3) {
                const neighbors = adj.get(i).sort((a, b) => atomVectors[i].distanceTo(atomVectors[a]) - atomVectors[i].distanceTo(atomVectors[b]));
                const toRemove = neighbors.slice(3);
                for(const neighbor of toRemove) {
                    bonds.delete(`${Math.min(i, neighbor)}-${Math.max(i, neighbor)}`);
                }
                changesMade = true;
            }
        }
        
        // Rebuild adjacency map after fixing over-coordinated atoms
        if (changesMade) {
            for(const [key, val] of adj) val.length = 0;
            for(const bondKey of bonds) {
                const [u,v] = bondKey.split('-').map(Number);
                adj.get(u).push(v);
                adj.get(v).push(u);
            }
        }

        for(let i=0; i<atomVectors.length; i++) {
            if(adj.get(i).length < 3) {
                const candidates = [];
                for(let j=0; j<atomVectors.length; j++) {
                   if(i === j || adj.get(i).includes(j)) continue;
                   candidates.push({index: j, distance: atomVectors[i].distanceTo(atomVectors[j]), coord: adj.get(j).length});
                }
                candidates.sort((a,b) => {
                   if (a.coord < 3 && b.coord >= 3) return -1;
                   if (a.coord >= 3 && b.coord < 3) return 1;
                   return a.distance - b.distance;
                });
                
                const addBond = (j, checkRing) => {
                    if (checkRing && checkIfFormsSmallRing(i, j, adj)) return false;
                    bonds.add(`${Math.min(i, j)}-${Math.max(i, j)}`);
                    adj.get(i).push(j);
                    adj.get(j).push(i);
                    changesMade = true;
                    return true;
                }

                // Pass 1
                for(const candidate of candidates) {
                    if(adj.get(i).length >= 3) break;
                    if(adj.get(candidate.index).length > 3) continue;
                    addBond(candidate.index, true);
                }
                // Pass 2
                if(adj.get(i).length < 3) {
                    for(const candidate of candidates) {
                        if(adj.get(i).length >= 3) break;
                        if(adj.get(candidate.index).length > 3) continue;
                        addBond(candidate.index, false); // Add bond even if it forms a small ring
                    }
                }
            }
        }
        if (!changesMade) break;
    }

    return Array.from(bonds).map(b => b.split('-').map(Number));
}
        
function checkIfFormsSmallRing(u, v, adj) {
    const uNeighbors = adj.get(u);
    const vNeighbors = adj.get(v);
    for(const n of uNeighbors) if(vNeighbors.includes(n)) return true;
    for (const un of uNeighbors) for (const vn of vNeighbors) if (adj.get(un).includes(vn)) return true;
    return false;
}

function relaxStructure(atomData, bonds, iterations, lengthStrength, angleStrength, pristineAngleStrength, planarityStrength, defectiveAtoms, idealLength, idealAngle) {
    let positions = atomData.map(a => new THREE.Vector3(...a.pos));
    const adj = new Map(positions.map((_, i) => [i, []]));
    for(const [u,v] of bonds) { adj.get(u).push(v); adj.get(v).push(u); }

    for (let i = 0; i < iterations; i++) {
        if (i % 20 === 0) postStatus(`Relaxing structure... (Iteration ${i + 1}/${iterations})`);
        const forces = positions.map(() => new THREE.Vector3(0, 0, 0));

        for (const [j, k] of bonds) {
            const p1 = positions[j], p2 = positions[k];
            const delta = new THREE.Vector3().subVectors(p2, p1);
            const distance = delta.length();
            if (distance > 1e-6) {
                const force = delta.normalize().multiplyScalar(lengthStrength * (distance - idealLength));
                forces[j].add(force);
                forces[k].sub(force);
            }
        }
        
        if (angleStrength > 0) {
            for (let j = 0; j < positions.length; j++) {
                const neighbors = adj.get(j);
                if (neighbors.length < 2) continue;
                for (let n1 = 0; n1 < neighbors.length; n1++) {
                    for (let n2 = n1 + 1; n2 < neighbors.length; n2++) {
                        const v1 = new THREE.Vector3().subVectors(positions[neighbors[n1]], positions[j]);
                        const v2 = new THREE.Vector3().subVectors(positions[neighbors[n2]], positions[j]);
                        if (v1.lengthSq() === 0 || v2.lengthSq() === 0) continue;
                        const angleDisplacement = v1.angleTo(v2) - idealAngle;
                        let finalForceMag = angleStrength * 10 * angleDisplacement * (defectiveAtoms.has(j) ? 1.0 : pristineAngleStrength);
                        const cross = new THREE.Vector3().crossVectors(v1, v2);
                        if (cross.lengthSq() > 1e-8) {
                            const fDir1 = new THREE.Vector3().crossVectors(v1, cross).normalize();
                            forces[neighbors[n1]].add(fDir1.multiplyScalar(finalForceMag));
                            const fDir2 = new THREE.Vector3().crossVectors(v2, cross.negate()).normalize();
                            forces[neighbors[n2]].add(fDir2.multiplyScalar(finalForceMag));
                        }
                    }
                }
            }
        }

        if (planarityStrength > 0) {
            for (let j = 0; j < positions.length; j++) {
                const neighbors = adj.get(j);
                if (neighbors.length !== 3) continue;
                const p_c = positions[j], p1 = positions[neighbors[0]], p2 = positions[neighbors[1]], p3 = positions[neighbors[2]];
                const normal = new THREE.Vector3().crossVectors(p2.clone().sub(p1), p3.clone().sub(p1));
                if (normal.lengthSq() < 1e-8) continue;
                normal.normalize();
                const dist = p_c.clone().sub(p1).dot(normal);
                const force = normal.clone().multiplyScalar(-planarityStrength * 5.0 * dist);
                forces[j].add(force);
                const reaction = force.clone().multiplyScalar(-1/3);
                forces[neighbors[0]].add(reaction);
                forces[neighbors[1]].add(reaction);
                forces[neighbors[2]].add(reaction);
            }
        }
        positions.forEach((pos, index) => pos.add(forces[index]));
    }
    return atomData.map((data, i) => ({ pos: [positions[i].x, positions[i].y, positions[i].z], tile: data.tile }));
}

function analyzeRings(bonds, numAtoms, edgeAtoms) {
    postStatus('Analyzing ring structure...');
    const adj = new Map(Array.from({length: numAtoms}, (_, i) => [i, []]));
    for (const [u, v] of bonds) { adj.get(u).push(v); adj.get(v).push(u); }

    const ringCounts = { 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, "9+":0 };
    const bondsInNonSixRings = new Set();
    const bfs = (start, end, excluded) => {
        const q = [[start, [start]]];
        const visited = new Set([start, ...excluded]);
        while (q.length > 0) {
            const [curr, path] = q.shift();
            if (curr === end) return path;
            for (const neighbor of adj.get(curr)) {
                if (!visited.has(neighbor)) {
                    visited.add(neighbor);
                    q.push([neighbor, [...path, neighbor]]);
                }
            }
        }
        return null;
    };

    for (const [u, v] of bonds) {
        if (edgeAtoms.has(u) || edgeAtoms.has(v)) continue;
        const adj_u = adj.get(u), adj_v = adj.get(v);
        adj.set(u, adj_u.filter(n => n !== v));
        adj.set(v, adj_v.filter(n => n !== u));
        let inNonSix = false;
        const path1 = bfs(u, v, new Set());
        if (path1) {
            if (path1.length !== 6) inNonSix = true;
            if (ringCounts[path1.length]) ringCounts[path1.length] += 1 / path1.length;
            else if (path1.length >= 9) ringCounts["9+"] += 1 / path1.length;
            const path2 = bfs(u, v, new Set(path1.slice(1,-1)));
            if(path2) {
                if (path2.length !== 6) inNonSix = true;
                if (ringCounts[path2.length]) ringCounts[path2.length] += 1 / path2.length;
                else if (path2.length >= 9) ringCounts["9+"] += 1 / path2.length;
            }
        } else {
            inNonSix = true;
        }
        if (inNonSix) bondsInNonSixRings.add(`${Math.min(u,v)}-${Math.max(u,v)}`);
        adj.set(u, adj_u);
        adj.set(v, adj_v);
    }
    return { bondsInNonSixRings: Array.from(bondsInNonSixRings), ringCounts };
}

const distToSegment2D = (p, segStart, segEnd) => {
    const l2 = segStart.distanceToSquared(segEnd);
    if (l2 === 0) return p.distanceTo(segStart);
    let t = ((p.x - segStart.x) * (segEnd.x - segStart.x) + (p.y - segStart.y) * (segEnd.y - segStart.y)) / l2;
    t = Math.max(0, Math.min(1, t));
    const proj = segStart.clone().lerp(segEnd, t);
    return p.distanceTo(proj);
};

function findEdgeAtoms(atomData) {
    postStatus('Identifying edge atoms...');
    const edgeAtoms = new Set();
    const edgeThreshold = 3.0;
    const atoms = atomData.map(data => new THREE.Vector3(...data.pos));
    if (atoms.length === 0) return edgeAtoms;

    const boundingBox = new THREE.Box3().setFromPoints(atoms);
    const center = boundingBox.getCenter(new THREE.Vector3());
    const furthestAtom = atoms.reduce((a, b) => center.distanceTo(a) > center.distanceTo(b) ? a : b);
    const modelRadius = center.distanceTo(furthestAtom);
    
    atoms.forEach((atom, i) => {
        if (center.distanceTo(atom) > modelRadius - edgeThreshold) {
            edgeAtoms.add(i);
        }
    });
    
    return edgeAtoms;
}


