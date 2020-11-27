import { setup_state, setup_opt } from "packing_wasm"

// Molecular buttons
const molRadius = document.getElementById("mol-radius")
const molDistance = document.getElementById("mol-distance")
const molAngle = document.getElementById("mol-angle")
const spaceGroup = document.getElementById("space-group")

// Interaction Buttons
const playPauseButton = document.getElementById("play-pause");
const resetButton = document.getElementById("reset-state");
const coolButton = document.getElementById("cool");

// Modifiable Values
const stepSize = document.getElementById("step-size");
const temperature = document.getElementById("temperature");
const iterationSteps = document.getElementById("iteration-steps");

// Output
const image = document.getElementById("output-figure");
const totalSteps = document.getElementById("total-steps");
const score = document.getElementById("score");

var cool;
var state;
var optimiser;

// Initialise state to default values
cool = false;
let animationId = null;

stepSize.value = 0.1
iterationSteps.value = 100n
molRadius.value = 0.7;
molDistance.value = 1.;
molAngle.value = 120.;
temperature.value = 0.2;
spaceGroup.value = "p2";
playPauseButton.textContent = "▶";

// Reset and reinitialise the state
function reset() {
    // These are value that need to be re-initialised when resetting the state
    totalSteps.value = 0n;
    temperature.value = 0.2;
    console.log(spaceGroup.value);
    state = setup_state(molRadius.value, molAngle.value, molDistance.value, spaceGroup.value);
    optimiser = setup_opt(temperature.value, stepSize.value, 1000n);
    image.innerHTML = state.as_svg();
    score.value = state.score().toFixed(3);
}

function renderLoop () {

    // Update values from the form
    optimiser.step_size = stepSize.value;
    optimiser.kt = Number(temperature.value);
    optimiser.steps = BigInt(iterationSteps.value);

    // Update and optimise state 
    state = optimiser.optimise_state(state);

    image.innerHTML = state.as_svg();
    stepSize.value = optimiser.step_size;
    totalSteps.value = BigInt(totalSteps.value) + optimiser.steps;
    score.value = state.score().toFixed(3);

    if (cool) {
        optimiser.kt *= 0.9;
    }
    if (stepSize.value < 1e-5) {
        pause();
        return
    }
    // console.log(optimiser.kt);

    temperature.value = optimiser.kt;

    // Request the next frame
    animationId = requestAnimationFrame(renderLoop);
};

const isPaused = () => {
    return animationId === null;
}

const play = () => {
  playPauseButton.textContent = "⏸";
  renderLoop();
};

const pause = () => {
  playPauseButton.textContent = "▶";
  cancelAnimationFrame(animationId);
  animationId = null;
};

playPauseButton.addEventListener("click", event => {
  if (isPaused()) {
    play();
  } else {
    pause();
  }
});

coolButton.addEventListener("click", event => {
  if (!cool) {
      cool = true
  }
});

resetButton.addEventListener("click", event => {
    reset();
});

reset();
// play();
