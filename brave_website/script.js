// JavaScript code (script.js)

// Function to toggle the visibility of the filter buttons
function toggleFilterPanel() {
    var filterPanel = document.getElementById('filter-panel');
    filterPanel.classList.toggle('show');
}

// Event listener to toggle the visibility of the filter buttons when "Apply Filters" button is clicked
document.getElementById('apply-filters-btn').addEventListener('click', function() {
    toggleFilterPanel();
});

// Function to toggle the visibility of the dropdown menu
function toggleDropdown(dropdownId) {
    var dropdown = document.getElementById(dropdownId);
    dropdown.classList.toggle('show');
}

// Event listener to handle clicks on filter buttons
document.addEventListener('DOMContentLoaded', function() {
    var filterButtons = document.querySelectorAll('.filter-btn');
    filterButtons.forEach(function(button) {
        button.addEventListener('click', function(event) {
            var dropdownId = event.target.dataset.dropdown;
            if (dropdownId) {
                toggleDropdown(dropdownId);
            }
        });
    });
});

// Close any open dropdown menus if the user clicks outside of them
window.addEventListener('click', function(event) {
    if (!event.target.matches('.filter-btn')) {
        var dropdowns = document.querySelectorAll('.dropdown-content');
        dropdowns.forEach(function(dropdown) {
            if (dropdown.classList.contains('show') && !dropdown.contains(event.target)) {
                dropdown.classList.remove('show');
            }
        });
    }
});

document.getElementById('search-form').addEventListener('submit', function(event) {
    event.preventDefault(); // Prevent form submission

    // Your search and filtering logic here...
});


